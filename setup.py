from setuptools import setup, find_packages


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="ecgem",
    version="0.0.1",
    description="Enzyme Constraint genome-scale models",
    keywords="python metabolic-models bioinformatics systems-biology",
    license="MIT",
    packages=find_packages(),
    install_requires=["cobra==0.20.0", "pandas==1.1.5"],
    zip_safe=False,
)
