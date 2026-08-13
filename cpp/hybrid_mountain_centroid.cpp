#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kMinimumHairpinLength = 3;
constexpr int kUnseen = -2;
constexpr int kUnpaired = -1;

bool CanPair(char left, char right) {
  return (left == 'A' && right == 'U') ||
         (left == 'U' && right == 'A') ||
         (left == 'G' && right == 'C') ||
         (left == 'C' && right == 'G') ||
         (left == 'G' && right == 'U') ||
         (left == 'U' && right == 'G');
}

std::string NormalizeSequence(std::string sequence) {
  if (sequence.empty()) {
    throw std::invalid_argument("Sequence must not be empty");
  }
  for (char& base : sequence) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    if (base == 'T') {
      base = 'U';
    }
    if (base != 'A' && base != 'C' && base != 'G' && base != 'U') {
      throw std::invalid_argument("Sequence contains an unsupported symbol");
    }
  }
  return sequence;
}

double MountainNormalizer(int length) {
  double total = 0.0;
  for (int cut = 1; cut < length; ++cut) {
    const double bound = std::min(cut, length - cut);
    total += bound * bound;
  }
  return total;
}

struct Layer {
  Layer(int sequence_length, int external_depth)
      : depth(external_depth),
        reduced_length(sequence_length - 2 * external_depth),
        costs(static_cast<std::size_t>(reduced_length) *
                  (reduced_length + 1) / 2,
              std::numeric_limits<double>::quiet_NaN()),
        decisions(costs.size(), kUnseen) {}

  std::size_t Index(int left, int right) const {
    const int row = left - depth;
    const int column = right - depth;
    if (row < 0 || column < row || column >= reduced_length) {
      throw std::logic_error("Infeasible interval-depth state");
    }
    const std::size_t prefix =
        static_cast<std::size_t>(row) * reduced_length -
        static_cast<std::size_t>(row) * (row - 1) / 2;
    return prefix + static_cast<std::size_t>(column - row);
  }

  int depth;
  int reduced_length;
  std::vector<double> costs;
  std::vector<int> decisions;
};

class HybridSolver {
 public:
  HybridSolver(const std::string& sequence,
               const std::vector<double>& expected_heights,
               const std::vector<double>& bpp,
               double alpha)
      : sequence_(sequence),
        expected_heights_(expected_heights),
        bpp_(bpp),
        alpha_(alpha),
        length_(static_cast<int>(sequence.size())),
        legal_partners_(length_),
        layers_((length_ + 1) / 2) {
    if (alpha_ < 0.0 || alpha_ > 1.0 || !std::isfinite(alpha_)) {
      throw std::invalid_argument("Alpha must lie in [0,1]");
    }
    if (expected_heights_.size() != static_cast<std::size_t>(length_ - 1) ||
        bpp_.size() != static_cast<std::size_t>(length_) * length_) {
      throw std::invalid_argument("Hybrid input dimensions do not match sequence");
    }
    const double profile_scale = MountainNormalizer(length_);
    const double pair_scale = length_ / 2;
    if (profile_scale <= 0.0 || pair_scale <= 0.0) {
      throw std::invalid_argument("Hybrid solver requires length at least two");
    }
    if (alpha_ < 1.0) {
      mountain_weight_ = 1.0;
      pair_weight_ = alpha_ * profile_scale /
                     ((1.0 - alpha_) * pair_scale);
    } else {
      mountain_weight_ = 0.0;
      pair_weight_ = 1.0;
    }
    for (int left = 0; left < length_; ++left) {
      for (int right = left + kMinimumHairpinLength + 1; right < length_;
           ++right) {
        if (CanPair(sequence_[left], sequence_[right])) {
          legal_partners_[left].push_back(right);
        }
      }
    }
  }

  void Run() {
    objective_ = Solve(0, length_ - 1, 0);
    structure_.assign(length_, '.');
    Traceback(0, length_ - 1, 0);
  }

  const std::string& structure() const { return structure_; }
  double objective() const { return objective_; }
  std::uint64_t states_evaluated() const { return states_evaluated_; }
  std::uint64_t partner_transitions_evaluated() const {
    return partner_transitions_evaluated_;
  }
  int effective_depth_levels() const { return maximum_external_depth_ + 1; }

 private:
  Layer& GetLayer(int depth) {
    if (depth < 0 || depth >= static_cast<int>(layers_.size())) {
      throw std::logic_error("External depth is out of range");
    }
    if (!layers_[depth]) {
      layers_[depth] = std::make_unique<Layer>(length_, depth);
    }
    return *layers_[depth];
  }

  double Bpp(int left, int right) const {
    return bpp_[static_cast<std::size_t>(left) * length_ + right];
  }

  double CutCost(int position, int depth) const {
    if (position == length_ - 1) {
      return 0.0;
    }
    const double difference = depth - expected_heights_[position];
    return mountain_weight_ * difference * difference;
  }

  double SubproblemCost(int left, int right, int depth) {
    return left > right ? 0.0 : Solve(left, right, depth);
  }

  double Solve(int left, int right, int depth) {
    if (depth > std::min(left, length_ - 1 - right)) {
      throw std::logic_error("Reached an infeasible external-depth state");
    }
    Layer& layer = GetLayer(depth);
    const std::size_t index = layer.Index(left, right);
    if (layer.decisions[index] != kUnseen) {
      return layer.costs[index];
    }

    ++states_evaluated_;
    maximum_external_depth_ = std::max(maximum_external_depth_, depth);
    double best_cost =
        CutCost(left, depth) + SubproblemCost(left + 1, right, depth);
    int best_partner = kUnpaired;

    for (int partner : legal_partners_[left]) {
      if (partner > right) {
        break;
      }
      ++partner_transitions_evaluated_;
      const double pair_gain = 2.0 * Bpp(left, partner) - 1.0;
      const double candidate =
          CutCost(left, depth + 1) +
          SubproblemCost(left + 1, partner - 1, depth + 1) +
          CutCost(partner, depth) +
          SubproblemCost(partner + 1, right, depth) -
          pair_weight_ * pair_gain;
      if (candidate < best_cost) {
        best_cost = candidate;
        best_partner = partner;
      }
    }

    layer.costs[index] = best_cost;
    layer.decisions[index] = best_partner;
    return best_cost;
  }

  void Traceback(int left, int right, int depth) {
    if (left > right) {
      return;
    }
    Layer& layer = GetLayer(depth);
    const int partner = layer.decisions[layer.Index(left, right)];
    if (partner == kUnpaired) {
      Traceback(left + 1, right, depth);
      return;
    }
    if (partner < 0) {
      throw std::logic_error("Missing traceback decision");
    }
    structure_[left] = '(';
    structure_[partner] = ')';
    Traceback(left + 1, partner - 1, depth + 1);
    Traceback(partner + 1, right, depth);
  }

  const std::string& sequence_;
  const std::vector<double>& expected_heights_;
  const std::vector<double>& bpp_;
  double alpha_;
  int length_;
  double mountain_weight_ = 0.0;
  double pair_weight_ = 0.0;
  std::vector<std::vector<int>> legal_partners_;
  std::vector<std::unique_ptr<Layer>> layers_;
  std::string structure_;
  double objective_ = 0.0;
  std::uint64_t states_evaluated_ = 0;
  std::uint64_t partner_transitions_evaluated_ = 0;
  int maximum_external_depth_ = 0;
};

}  // namespace

int main() {
  try {
    std::string sequence;
    if (!std::getline(std::cin, sequence)) {
      throw std::invalid_argument("Missing sequence line");
    }
    sequence = NormalizeSequence(std::move(sequence));
    const int length = static_cast<int>(sequence.size());

    int alpha_count = 0;
    if (!(std::cin >> alpha_count) || alpha_count < 1) {
      throw std::invalid_argument("Missing or invalid alpha count");
    }
    std::vector<double> alphas(alpha_count);
    for (double& alpha : alphas) {
      if (!(std::cin >> alpha)) {
        throw std::invalid_argument("Missing alpha value");
      }
    }

    std::vector<double> expected_heights(std::max(0, length - 1));
    for (double& value : expected_heights) {
      if (!(std::cin >> value) || !std::isfinite(value)) {
        throw std::invalid_argument("Missing or invalid expected mountain height");
      }
    }

    std::vector<double> bpp(static_cast<std::size_t>(length) * length, 0.0);
    for (int left = 0; left < length; ++left) {
      for (int right = left + 1; right < length; ++right) {
        double probability = 0.0;
        if (!(std::cin >> probability) || probability < 0.0 ||
            probability > 1.0 || !std::isfinite(probability)) {
          throw std::invalid_argument("Missing or invalid BPP value");
        }
        bpp[static_cast<std::size_t>(left) * length + right] = probability;
      }
    }

    std::cout << std::setprecision(17);
    for (double alpha : alphas) {
      const auto started = std::chrono::steady_clock::now();
      HybridSolver solver(sequence, expected_heights, bpp, alpha);
      solver.Run();
      const std::chrono::duration<double> elapsed =
          std::chrono::steady_clock::now() - started;
      std::cout << alpha << '\t' << solver.structure() << '\t'
                << solver.objective() << '\t' << solver.states_evaluated()
                << '\t' << solver.partner_transitions_evaluated() << '\t'
                << solver.effective_depth_levels() << '\t' << elapsed.count()
                << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "hybrid_mountain_centroid: " << error.what() << '\n';
    return 2;
  }
}
