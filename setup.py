from setuptools import setup, find_packages

setup(
    name="julia2",
    version="0.1.0",
    description="A command-line tool for detecting index hopping in RNA sequence data",
    author="Benjamin Glick",
    author_email="glick@glick.cloud",  # Replace with your email
    url="https://github.com/benhg/julia2-tool",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # Add your dependencies here, e.g., "numpy", "pandas", etc.
    ],
    entry_points={
        "console_scripts": [
            "julia2=julia2:main",  # Maps "julia2" executable to the "main" function in julia2.py
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
