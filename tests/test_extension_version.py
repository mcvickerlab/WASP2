import importlib.metadata

import wasp2_rust


def test_rust_extension_version_matches_distribution() -> None:
    assert wasp2_rust.__version__ == importlib.metadata.version("wasp2")
