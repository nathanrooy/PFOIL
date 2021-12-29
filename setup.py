from numpy.distutils.core import Extension

fortran_ext = Extension(name = 'xfoil', sources = ['pfoil/xfoil.f'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(
        name         = 'pfoil',
        version      = '0.0.3',
        author       = "Nathan A. Rooy",
        author_email = "nathanrooy@gmail.com",
        description  = "XFOIL based aerodynamic analysis and design",
        ext_modules  = [fortran_ext],
        packages     = ['pfoil'],
        install_requires = [
            'numpy>=1.19'
        ]
    )
