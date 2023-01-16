from setuptools import setup

setup(
    install_requires=['pandas','argparse','matplotlib'],
    entry_points={
        'console_scripts': [
            'tRNAmap = tRNAmap:main'
        ]
    }
)

#pip install -e .

