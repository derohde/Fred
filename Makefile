install:
	python setup.py install --user

clean:
	rm -r dist build/ Fred_Frechet.egg-info/ & pip uninstall Fred-Frechet -y

