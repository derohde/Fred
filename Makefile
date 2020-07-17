pre:
	sudo apt install -y libboost-all-dev
	sudo apt-get install -y python3-setuptools
	sudo apt-get install -y python3-numpy
	sudo apt-get install -y python3-matplotlib
	sudo apt-get install -y cmake

install:
	cd py && /usr/bin/python3 ./setup.py install --user

clean:
	rm -r py/dist py/build/ py/Fred.egg-info/
	pip3 uninstall Fred -y

