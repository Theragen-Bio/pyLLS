from setuptools import setup, find_packages


setup(
    name='pyLLS',
    version='0.5',
    license='MIT',
    author="Sejin Oh",
    author_email='agicic@naver.com',
    packages=find_packages('pyLLS'),
    #package_dir={'': 'pyLLS'},
    url='https://github.com/osj118/pyLLS',
    keywords='missing value imputation',
    install_requires=[
          'scikit-learn',
          'pandas',
          'numpy',
          'kneed',
          'tqdm'
      ],
    python_requires='>=3.7',

)