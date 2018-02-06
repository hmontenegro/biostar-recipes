USER=www
DATA_DIR=/export/sites/main_data/initial
DATA_HOST=data.bioinformatics.recipes

# The initial data for all recipes
DATA_FILE=recipes-initial-data.tar.gz

all:
	@echo Valid commands: push, data_pull

# The local directory
export/local:
	mkdir -p export/local

# Make the required directories
dir: export/local

push:
	git commit -am "Update by `whoami` on `date` from `hostname`"
	git push

mothur:
	python manage.py project --json projects/metagenome/mothur-project.hjson --privacy public --jobs

data: dir
	(cd export && curl http://data.bioinformatics.recipes/initial/${DATA_FILE} > ${DATA_FILE} )
	(cd export && tar zxvf ${DATA_FILE})

pack: dir
	(cd export && tar czvf ${DATA_FILE} local )
	(cd export && rsync -avz ${DATA_FILE} ${USER}@${DATA_HOST}:${DATA_DIR}/)

fish:
	python manage.py project --json initial/fish-project.hjson

giraffe:
	python manage.py project --root ../biostar-recipes --json projects/giraffe/giraffe-project.hjson --sticky --privacy public

mothur:
	python manage.py project --root ../biostar-recipes --json projects/metagenome/mothur-project.hjson --privacy public
