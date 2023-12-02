import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="aberration_multigraph",
    version="0.0.1",
    author="Siddharth Sheth",
    author_email="sheth.sid@gmail.com",
    description="A package to model aberration multigraphs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/siddharthsheth/aberration_multigraph",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=['networkx'],
    # extras_require={
    #     'dev': [
    #             'sphinx',
    #             'sphinx-rtd-theme',
    #             'sphinxcontrib-bibtex',
    #             'sphinx-click',
    #             'pytest',
    #             'pycodestyle',
    #             'autopep8',
    #             'pytest-cov',
    #            ],
    #     'viz': [
    #             'ds2viz'
    #            ]
    #     },
    # entry_points={
    #     'console_scripts': [
    #         'greedypermutation=greedypermutation.cli:cli',
    #     ],
    # },
)