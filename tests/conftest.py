import os
import pytest
from pathlib import Path

@pytest.fixture(scope="module")
def get_current_dir():
    return str(os.path.abspath(os.path.join(os.path.dirname(__file__))))

@pytest.fixture(scope="module")
def get_parent_dir(get_current_dir):
    return Path(get_current_dir).parent

@pytest.fixture(scope="module")
def get_steinbock_out_dir_mcd(get_parent_dir):
    return os.path.join(get_parent_dir, 'data', 'test_mcd')

@pytest.fixture(scope="module")
def get_steinbock_out_dir_txt(get_parent_dir):
    return os.path.join(get_parent_dir, 'data', 'test_txt')

@pytest.fixture(scope="module")
def get_steinbock_out_dir_tiff(get_parent_dir):
    return os.path.join(get_parent_dir, 'data', 'test_tiff')