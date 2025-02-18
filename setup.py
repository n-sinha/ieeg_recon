from setuptools import setup, find_packages

# Read requirements
def read_requirements(filename):
    """Read requirements from requirements.txt"""
    with open(filename) as f:
        requirements = []
        extras_require = {'dev': []}
        
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            if '; extra == ' in line:
                # This is a development requirement
                req = line.split('; extra == ')[0].strip()
                extras_require['dev'].append(req)
            else:
                # This is a core requirement
                requirements.append(line)
                
        return requirements, extras_require

requirements, extras_require = read_requirements('requirements.txt')

setup(
    # Basic package information
    name="ieeg_recon",
    version="1.0.0",
    author="Nishant Sinha",
    author_email="nishants@seas.upenn.edu",
    description="A pipeline for reconstructing iEEG electrode locations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    
    # Project URLs
    url="https://github.com/n-sinha/ieeg_recon",
    
    # Find packages automatically (looks for __init__.py files)
    packages=find_packages(include=['ieeg_recon', 'ieeg_recon.*']),
    
    # Package dependencies
    install_requires=requirements,
    
    # Development dependencies (optional)
    extras_require=extras_require,
    
    # Python version compatibility
    python_requires='>=3.8',
    
    # Package classifiers (helps people find your package)
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
    ],
    
    # Include non-Python files
    package_data={
        'src': ['*.txt', '*.json'],
    },
    
    # Command-line scripts
    entry_points={
        'console_scripts': [
            'ieeg-recon=src.run_ieeg_recon:main',
        ],
    },
    
    # Test suite
    test_suite='tests',
) 