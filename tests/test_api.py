from types import SimpleNamespace

import numpy as np

import mountain_centroid.api as api
from mountain_centroid.api import predict_from_profile, predict_hybrid
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


def test_public_prediction_defaults_to_vienna(monkeypatch):
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

    assert captured["backend"] == "vienna"
    assert prediction.bpp_backend == "vienna"


def test_public_hybrid_curve_computes_bpp_once(monkeypatch):
    captured = {"bpp_calls": 0}
    bpp = np.zeros((5, 5), dtype=float)
    mu = [0.0] * 4

    def fake_compute_bpp_and_mu(sequence, **kwargs):
        captured["bpp_calls"] += 1
        captured["backend"] = kwargs["backend"]
        return bpp, mu

    def fake_cpp_hybrid(sequence, expected_heights, probabilities, alphas, **kwargs):
        captured["alphas"] = tuple(alphas)
        assert expected_heights is mu
        assert probabilities is bpp
        return tuple(
            SimpleNamespace(
                alpha=alpha,
                structure=".....",
                pairs=(),
                heights=(0, 0, 0, 0, 0, 0),
                squared_mountain_error=0.0,
                centroid_gain=0.0,
                hybrid_objective=0.0,
            )
            for alpha in alphas
        )

    monkeypatch.setattr(api, "compute_bpp_and_mu", fake_compute_bpp_and_mu)
    monkeypatch.setattr(api, "cpp_hybrid_mountain_centroid", fake_cpp_hybrid)
    predictions = api.predict_hybrid_curve("GAAAU", (0.0, 0.25, 1.0))

    assert captured["bpp_calls"] == 1
    assert captured["backend"] == "vienna"
    assert captured["alphas"] == (0.0, 0.25, 1.0)
    assert tuple(prediction.alpha for prediction in predictions) == captured["alphas"]
    assert all(prediction.bpp_backend == "vienna" for prediction in predictions)


def test_public_single_hybrid_prediction(monkeypatch):
    expected = object()

    def fake_curve(sequence, alphas, **kwargs):
        assert sequence == "GAAAU"
        assert alphas == (0.2,)
        return (expected,)

    monkeypatch.setattr(api, "predict_hybrid_curve", fake_curve)
    assert predict_hybrid("GAAAU", 0.2) is expected
