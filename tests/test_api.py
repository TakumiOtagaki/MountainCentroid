import mountain_centroid.api as api
from mountain_centroid.api import predict_from_profile
from mountain_centroid.constrained import sequence_constrained_mountain_centroid


def _python_solver(sequence, expected_heights, *, executable=None):
    return sequence_constrained_mountain_centroid(sequence, expected_heights)


def test_public_profile_api_normalizes_sequence_and_reports_direct_error(
    monkeypatch,
):
    monkeypatch.setattr(
        api,
        "cpp_sequence_constrained_mountain_centroid",
        _python_solver,
    )
    mu = [1.0, 1.0, 1.0, 1.0]
    prediction = predict_from_profile("gaaAt", mu)

    direct_error = sum(
        (height - expected) ** 2
        for height, expected in zip(
            prediction.mountain_heights[1:-1],
            prediction.expected_mountain_heights,
        )
    )
    assert prediction.sequence == "GAAAU"
    assert prediction.structure == "(...)"
    assert prediction.bpp_backend == "provided"
    assert prediction.solver_backend == "cpp_interval_depth_dp"
    assert prediction.squared_mountain_error == direct_error


def test_public_prediction_defaults_to_linearpartition(monkeypatch):
    captured = {}

    def fake_compute_bpp_and_mu(sequence, **kwargs):
        captured.update(kwargs)
        return None, [0.0] * (len(sequence) - 1)

    monkeypatch.setattr(api, "compute_bpp_and_mu", fake_compute_bpp_and_mu)
    monkeypatch.setattr(
        api,
        "cpp_sequence_constrained_mountain_centroid",
        _python_solver,
    )
    prediction = api.predict("GAAAU")

    assert captured["backend"] == "linearpartition"
    assert prediction.bpp_backend == "linearpartition"
