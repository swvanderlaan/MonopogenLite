import setuptools
from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MonopogenLite",
    version="v1.0.0",
    author="Jinzhuang Dou; Sander W. van der Laan",
    author_email="jdou1@mdanderson.org; s.w.vanderlaan-2@umcutrecht.nl",
    description="SNV calling from single cell sequencing",
    long_description=long_description,
    long_description_content_type="",
    url="",
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL v3 or later",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)

