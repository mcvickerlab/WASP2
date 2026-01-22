"""Tests for the WASP2 CLI module."""

import pytest
from typer.testing import CliRunner
from wasp2.cli import app
from wasp2.config import WASP2Config, load_config, get_config_path


runner = CliRunner()


class TestCliBasics:
    """Test basic CLI functionality."""

    def test_help_shows_commands(self):
        """Test that --help shows all main commands."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "count" in result.output
        assert "analyze" in result.output
        assert "map" in result.output
        assert "info" in result.output
        assert "config" in result.output

    def test_version_flag(self):
        """Test --version flag displays version."""
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "WASP2" in result.output
        assert "version" in result.output.lower()

    def test_verbosity_flag_accepted(self):
        """Test -v flag is accepted."""
        result = runner.invoke(app, ["-v", "--help"])
        assert result.exit_code == 0


class TestInfoCommand:
    """Test the info command."""

    def test_info_shows_system_info(self):
        """Test info command shows system information."""
        result = runner.invoke(app, ["info"])
        assert result.exit_code == 0
        assert "WASP2 Version" in result.output
        assert "Python Version" in result.output
        assert "Platform" in result.output

    def test_info_shows_dependencies(self):
        """Test info command shows dependency status."""
        result = runner.invoke(app, ["info"])
        assert result.exit_code == 0
        assert "Dependencies" in result.output
        assert "numpy" in result.output
        assert "pandas" in result.output

    def test_info_shows_config(self):
        """Test info command shows configuration."""
        result = runner.invoke(app, ["info"])
        assert result.exit_code == 0
        assert "Configuration" in result.output
        assert "Config File" in result.output


class TestConfigCommand:
    """Test the config command."""

    def test_config_show(self):
        """Test config show displays settings."""
        result = runner.invoke(app, ["config", "show"])
        assert result.exit_code == 0
        assert "threads" in result.output
        assert "use_rust" in result.output

    def test_config_path(self):
        """Test config path shows file location."""
        result = runner.invoke(app, ["config", "path"])
        assert result.exit_code == 0
        assert ".wasp2" in result.output
        assert "config.yaml" in result.output

    def test_config_help(self):
        """Test config subcommand help."""
        result = runner.invoke(app, ["config", "--help"])
        assert result.exit_code == 0
        assert "show" in result.output
        assert "set" in result.output
        assert "reset" in result.output


class TestCountCommand:
    """Test the count subcommand."""

    def test_count_help(self):
        """Test count --help shows commands."""
        result = runner.invoke(app, ["count", "--help"])
        assert result.exit_code == 0
        assert "count-variants" in result.output
        assert "count-variants-sc" in result.output

    def test_count_variants_help(self):
        """Test count-variants --help shows options."""
        result = runner.invoke(app, ["count", "count-variants", "--help"])
        assert result.exit_code == 0
        assert "BAM" in result.output
        assert "--dry-run" in result.output
        assert "--samples" in result.output


class TestMapCommand:
    """Test the map subcommand."""

    def test_map_help(self):
        """Test map --help shows commands."""
        result = runner.invoke(app, ["map", "--help"])
        assert result.exit_code == 0
        assert "make-reads" in result.output
        assert "filter-remapped" in result.output

    def test_make_reads_help(self):
        """Test make-reads --help shows options."""
        result = runner.invoke(app, ["map", "make-reads", "--help"])
        assert result.exit_code == 0
        assert "BAM" in result.output
        assert "--dry-run" in result.output


class TestAnalyzeCommand:
    """Test the analyze subcommand."""

    def test_analyze_help(self):
        """Test analyze --help shows commands."""
        result = runner.invoke(app, ["analyze", "--help"])
        assert result.exit_code == 0
        assert "find-imbalance" in result.output


class TestConfigModule:
    """Test the config module directly."""

    def test_default_config_values(self):
        """Test default configuration values."""
        config = WASP2Config()
        assert config.threads == 1
        assert config.use_rust is True
        assert config.output_format == "tsv"
        assert config.log_level == "WARNING"

    def test_config_to_dict(self):
        """Test config to_dict method."""
        config = WASP2Config()
        d = config.to_dict()
        assert "threads" in d
        assert "use_rust" in d
        assert d["threads"] == 1

    def test_get_config_path_returns_path(self):
        """Test get_config_path returns a Path."""
        path = get_config_path()
        assert path.suffix == ".yaml"
        assert ".wasp2" in str(path)
