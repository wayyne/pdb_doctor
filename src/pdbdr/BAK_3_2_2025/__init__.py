# pdbdr/__init__.py

"""
pdbdr: A package that implements fill_and_transfer with local or Forge ESM3.
No logic is changed from the monolithic version, only reorganized into modules.
"""

__version__ = "0.1.0"

from .pdb_io import *
from .esm_filling import *
from .alignment import *
from .stitching import *
from .sanity_check import *
from .constants import *
from tqdm import tqdm
import warnings

# Store the original __init__ function before overriding
if not hasattr(tqdm, "_original_init"):
    tqdm._original_init = tqdm.__init__

# Define default tqdm settings globally
tqdm_defaults = {
    "ncols": 80,  # Fixed width
    "ascii": True,  # Use ASCII progress bars
    "leave": False   # Keep progress bar after completion
}

# Override tqdm globally with safe recursion handling
def custom_tqdm_init(self, *args, **kwargs):
    kwargs = {**tqdm_defaults, **kwargs}
    self._original_init(*args, **kwargs)

tqdm.__init__ = custom_tqdm_init

# Suppress warnings about missing metadata; we know what we're doing!
warnings.filterwarnings(
    "ignore",
    message="Entity ID not found in metadata, using None as default"
)

warnings.filterwarnings("ignore", category=FutureWarning, module="esm.models.vqvae")
