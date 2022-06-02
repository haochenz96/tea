import setuptools 
import pathlib
import pkg_resources

# fetch from requirements.txt
with pathlib.Path("requirements.txt").open() as requirements_txt:
    install_requires = [
        str(requirement)
        for requirement in pkg_resources.parse_requirements(requirements_txt)
    ]

setuptools.setup(
    name='tea',
    version='0.1.0',
    author = 'Haochen Zhang',
    author_email = 'zhangh5@mskcc.org',
    description = 'a set of tools to do Tapestri Express Analysis',
    long_description = 'README.md',
    long_description_content_type = 'text/markdown',
    project_urls = {
        'Source': 'https://github.com/haochenz96/tea',
        'Bug Tracker': 'https://github.com/haochenz96/tea/issues',
    },
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    package_dir={"": "src"},
    packages = setuptools.find_packages(where="src"),
    install_requires=install_requires,
    # entry_points={
    #     'console_scripts': [
    #         'tapestri = cli.__main__:main',
    #     ],
    # },
)