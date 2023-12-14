from setuptools import setup

with open("README", 'r') as f:
    long_description = f.read()

setup(
   name='foo',
   version='1.0',
   description='A useful module',
   license="MIT",
   long_description=long_description,
   author='Man Foo',
   author_email='foomail@foo.example',
   url="http://www.foopackage.example/",
   packages=['foo'],  #same as name
   install_requires=['numpy', 'matplotlib', 'scipy', 'gmsh', 'pandas'], #external packages as dependencies
   scripts=[
            'scripts/cool',
            'scripts/skype',
           ]
)
