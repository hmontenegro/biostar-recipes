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

data: dir
	(cd export && curl http://${DATA_HOST}/initial/${DATA_FILE} > ${DATA_FILE} )
	(cd export && tar zxvf ${DATA_FILE})

pack: dir
	(cd export && tar czvf ${DATA_FILE} local )
	(cd export && rsync -avz ${DATA_FILE} ${USER}@${DATA_HOST}:${DATA_DIR}/)
