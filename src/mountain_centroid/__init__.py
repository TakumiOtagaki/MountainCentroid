"""Mountain Centroid public package API."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mountain-centroid")
except PackageNotFoundError:
    __version__ = "0+unknown"

from .api import (
    HybridPrediction,
    Prediction,
    predict,
    predict_from_profile,
    predict_hybrid,
    predict_hybrid_curve,
)
from .evaluation import (
    base_pair_f1,
    mean_squared_mountain_distance,
    normalized_squared_mountain_distance,
    squared_mountain_distance,
)
from .constrained import (
    ConstrainedDiagnostics,
    ConstrainedResult,
    sequence_constrained_mountain_centroid,
)
from .cpp_constrained import cpp_sequence_constrained_mountain_centroid
from .cpp_hybrid import HybridResult, cpp_hybrid_mountain_centroid
from .relaxed import RelaxedResult, relaxed_mountain_centroid

__all__ = [
    "__version__",
    "Prediction",
    "HybridPrediction",
    "HybridResult",
    "ConstrainedDiagnostics",
    "ConstrainedResult",
    "RelaxedResult",
    "base_pair_f1",
    "cpp_sequence_constrained_mountain_centroid",
    "cpp_hybrid_mountain_centroid",
    "mean_squared_mountain_distance",
    "normalized_squared_mountain_distance",
    "predict",
    "predict_from_profile",
    "predict_hybrid",
    "predict_hybrid_curve",
    "relaxed_mountain_centroid",
    "sequence_constrained_mountain_centroid",
    "squared_mountain_distance",
]
