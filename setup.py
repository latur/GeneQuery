from setuptools import setup, Extension

setup(
    name = 'genequery',
    version = '1.5.1',
    description = 'Expression-based phenotype searching engine',
    url = "https://github.com/latur/GeneQuery",
    author = 'Igor V.',
    author_email = 'latur@me.com',
    ext_modules = [Extension(
      'genequery', ['src/genequery.c'],
      include_dirs = ['src'],
      extra_compile_args=['-std=c99', '-m64', '-O3'] # '-lpthread'
    )],
    py_modules = ['genequery']
)

# rm -rf genequery.* dist build
# python3 setup.py build
# python3 setup.py install
# rm -rf genequery.* dist build && python3 setup.py build && python3 setup.py install
