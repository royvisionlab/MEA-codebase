from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

import sys

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    flags = ['-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')



# only bother to deal with compilation on Unix and MacOS
def make_pybind11_extension_with_flags (module_name, dependencies):

    c_opts = ['-O3', '-ffast-math']
    l_opts = []

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-std=c++11']

        c_opts = c_opts + darwin_opts
        l_opts = l_opts + darwin_opts

    if sys.platform == 'linux':
        linux_additional_opts = ['-fopenmp', ]

        c_opts = c_opts + linux_additional_opts
        l_opts = l_opts + linux_additional_opts

    return Extension(
        module_name,
        dependencies,
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True)
        ],
        language='c++',
        extra_compile_args = c_opts,
        extra_link_args= l_opts
    )


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': ['-O3', '-ffast-math'],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        print(opts)

        build_ext.build_extensions(self)

__version__ = '0.2'

py_modules = [
        'analysis_utils',
        'cell_typing',
        'config',
        'correlations',
        'lnp',
        'symphony_data',
        'vision_utils'
]

pybind11_ext_modules = [
    make_pybind11_extension_with_flags('correlations.correlation_cpp_extensions',
                                       ['correlations/cpp_extensions/ccg.cpp']),
    make_pybind11_extension_with_flags('vision_utils.vision_util_cpp_extensions',
                                       ['vision_utils/cpp_extensions/vision_util_extensions.cpp']),
]

setup(
    name='analysis',
    version=__version__,
    description='Utility functions for analysing MEA data collected in Symphony-DAS',
    author='Mike Manookin',
    # scripts=['config.py','lnp.py','sorting_hat.py','symphony_analysis.py','symphony_data.py','vision_utils.py'],
    # packages=find_packages(),
    install_requires=['pybind11>=2.3'],
    setup_requires=['pybind11>=2.3'],
    packages=py_modules,
    ext_modules=pybind11_ext_modules
    #cmdclass={'build_ext_cpp': BuildExt},
)
