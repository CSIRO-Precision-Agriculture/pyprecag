SHELL=/bin/sh
PACKAGE_NAME=pyprecag

.PHONY: help
help:
	echo
	echo 'Utility Makefile for pyprecag'
	echo '============================'
	echo
	echo 'Targets supported are:'
	echo
	echo '  * clean: removes the build directories, as well as __pycache__ and *.pyc files. Note that a clean also removes the generated documentation (as this is placed into build/docs).'
	echo '  * sitepkg-develop: Use this only when you are using a virtual environment created with the --system-site-packages flag. It forcibly installs some packages into the virtual environment to work around issues where a package and its plugins are installed to different locations'
	echo '  * install-deps: installs development and test dependencies into your virtual environment.'
	echo '  * develop: installs pyprecag in development mode.'
	echo '  * uninstall: removes the development package from pip.'
	echo '  * test: runs all unit tests.'
	echo '  * twine_check: Run `twine check`.'
	echo '  * lint: runs pylint.'
	echo '  * html: builds the HTML documentation.'
	echo '  * pdf: builds the documentation in PDF format.'
	echo '  * latex: builds LaTeX source, used to generate other formats.'
	echo '  * alldocs: builds all documentation formats.'
	echo '  * sdist: builds a source distribution.'
	echo '  * bdist_wheel: builds a universal wheel distribution.'
	echo '  * upload: uploads the source distribution and wheels to PyPI'
	echo '  * upload-test: uploads a test version to https://test.pypi.org/project/pyprecag'

.PHONY: twine_check
twine_check: sdist
	twine check ./dist/*

.PHONY: test
test:
	python -m unittest discover -s $(PACKAGE_NAME)/tests -p "*test*.py" -v

.PHONY: clean
clean:
	echo Cleaning ...
	rm -rf build/
	-find ./$(PACKAGE_NAME)/ -name "__pycache__" -exec rm -rf {} \;
	-find ./$(PACKAGE_NAME)/ -name "*.pyc" -exec rm -rf {} \;
	echo ... done

.PHONY: install-deps
install-deps:
	pip install -e.[dev,test]

.PHONY: develop
develop: install-deps
	python setup.py develop

.PHONY: uninstall
uninstall:
	-pip uninstall --yes $(PACKAGE_NAME)
	rm -rf *.egg-info/

#--system-site-packages plugin issue workaround
# To work correctly with virtual environments, sphinx, pytest, and all their
# plugins must be installed into the virtual environment even if they are
# already available in the site packages.
# If additional plugins to sphinx and pytest are added to setup.py, they must
# also be added here.
.PHONY: sitepkg-develop
sitepkg-develop: develop
	pip install --ignore-installed \
	sphinx \
	pytest pytest-cov pytest-sugar \
	ipython ipdb

.PHONY: lint
lint:
	pylint ./$(PACKAGE_NAME)/

.PHONY: sdist
sdist:
	python setup.py sdist

.PHONY: bdist_wheel
bdist_wheel:
	python setup.py bdist_wheel

.PHONY: html
html:
	sphinx-build -b html docs build/docs/html

.PHONY: latex
latex:
	sphinx-build -b latex docs build/docs/latex

pdf: latex
	$(MAKE) -C build/docs/latex all-pdf
	mkdir -p ./build/docs/pdf/
	mv build/docs/latex/*.pdf build/docs/pdf/

alldocs: html latex pdf

.PHONY: upload
upload: clean
	rm -rf dist/
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: upload-test
upload-test: clean
		rm -rf dist/
		python setup.py sdist bdist_wheel
		twine upload --repository-url https://test.pypi.org/legacy/ dist/*
