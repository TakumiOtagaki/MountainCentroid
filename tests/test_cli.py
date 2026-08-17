from types import SimpleNamespace

import pytest

from mountain_centroid import __version__
import mountain_centroid.mountain_pipeline as cli


def test_cli_reports_package_version(monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["mountain-centroid", "--version"])

    with pytest.raises(SystemExit) as exit_info:
        cli.main()

    assert exit_info.value.code == 0
    assert capsys.readouterr().out == f"mountain-centroid {__version__}\n"


def test_cli_without_alpha_uses_mountain_centroid(monkeypatch, capsys):
    prediction = SimpleNamespace(
        sequence="GAAAU",
        structure=".....",
        expected_mountain_heights=(0.0, 0.0, 0.0, 0.0),
        squared_mountain_error=0.0,
        bpp_backend="vienna",
    )
    monkeypatch.setattr(cli, "predict", lambda sequence, **kwargs: prediction)
    monkeypatch.setattr("sys.argv", ["mountain-centroid", "--seq", "GAAAU"])

    cli.main()

    assert "dot_bracket = ....." in capsys.readouterr().out


def test_cli_with_alphas_uses_one_hybrid_curve_call(monkeypatch, capsys):
    captured = {}

    def fake_curve(sequence, alphas, **kwargs):
        captured["sequence"] = sequence
        captured["alphas"] = alphas
        return tuple(
            SimpleNamespace(
                sequence="GAAAU",
                alpha=alpha,
                structure=".....",
                squared_mountain_error=0.0,
                centroid_gain=0.0,
                hybrid_objective=0.0,
                bpp_backend="vienna",
            )
            for alpha in alphas
        )

    monkeypatch.setattr(cli, "predict_hybrid_curve", fake_curve)
    monkeypatch.setattr(
        "sys.argv",
        [
            "mountain-centroid",
            "--seq",
            "GAAAU",
            "--alpha",
            "0",
            "--alpha",
            "0.23",
            "--alpha",
            "1",
        ],
    )

    cli.main()

    output = capsys.readouterr().out
    assert captured == {"sequence": "GAAAU", "alphas": [0.0, 0.23, 1.0]}
    assert output.count("dot_bracket = .....") == 3
    assert "alpha = 0.23" in output
