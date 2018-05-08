
from setuptools import setup
from setuptools_rust import RustExtension, Binding

setup(name='rust_fastforward',
      version='0.0.1',
      rust_extensions=[
          RustExtension('rust_fastforward._rust_fastforward', 'extensions/Cargo.toml',
                        binding=Binding.RustCPython)],
      packages=['rust_fastforward'],
      zip_safe=False)
