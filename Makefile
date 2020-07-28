# Data setup

data/hsal:
	cd data/raw/hsal; python parse.py

data/paer:
	cd data/raw/paer; python parse.py

data: data/paer data/hsal
.PHONY: data data/paer data/hsal

# Posterior inference

models := mnull mbatch mfull

# core recipes
sampling/hsalinarum/combined/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p hsalinarum/combined/$*/;\
	python hsalinarum.py $(word 4, $(subst /, ,$(@D))) $(word 5, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=20 \
	> hsalinarum/combined/$*/log.out \
	2> hsalinarum/combined/$*/log.err

sampling/hsalinarum/individual/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p hsalinarum/individual/$*/;\
	python hsalinarum.py $(word 4, $(subst /, ,$(@D))) mnull \
	--dataset $(word 5, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=20 \
	> hsalinarum/individual/$*/log.out \
	2> hsalinarum/individual/$*/log.err

sampling/paeruginosa/combined/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p paeruginosa/combined/$*/;\
	python paeruginosa.py $(word 4, $(subst /, ,$(@D))) $(word 5, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=15 \
	> paeruginosa/combined/$*/log.out \
	2> paeruginosa/combined/$*/log.err

sampling/paeruginosa/individual/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p paeruginosa/individual/$*/;\
	python paeruginosa.py $(word 4, $(subst /, ,$(@D))) mnull \
	--dataset $(word 5, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=15 \
	> paeruginosa/individual/$*/log.out \
	2> paeruginosa/individual/$*/log.err

# H. salinarum
conditions := low hi
hbatches := $(foreach batch,$(wildcard data/hi-oxidative/*),$(lastword $(subst /, ,$(batch))))

hsalCombined: $(foreach cond,$(conditions),$(foreach model,$(models),sampling/hsalinarum/combined/$(cond)/$(model)/samples/posterior_0.pkl))
hsalIndiv: $(foreach cond,$(conditions),$(foreach batch,$(hbatches),sampling/hsalinarum/individual/$(cond)/$(batch)/samples/posterior_0.pkl))
hsal: hsalCombined hsalIndiv
.PHONY: hsalCombined hsalIndiv hsal

# P. aeruginosa
acids := benzoate citric malic

paerCombined: $(foreach acid,$(acids),$(foreach model,$(models),sampling/paeruginosa/combined/$(acid)/$(model)/samples/posterior_0.pkl))
paerIndiv: sampling/paeruginosa/individual/benzoate/PA01_Benzoate_15_min_time_points/samples/posterior_0.pkl sampling/paeruginosa/individual/benzoate/PA01_Benzoate_repeat_19.07.17/samples/posterior_0.pkl sampling/paeruginosa/individual/malic/PA01_Malic_09.03.17/samples/posterior_0.pkl sampling/paeruginosa/individual/malic/PA01_Malic_repeat_27.07.17/samples/posterior_0.pkl sampling/paeruginosa/individual/citric/PA01_citric_15_min_time_points_06.03.17/samples/posterior_0.pkl sampling/paeruginosa/individual/citric/PA01_Citric_rerun_11.07.17/samples/posterior_0.pkl

paer: paerCombined paerIndiv
.PHONY: paer paerCombined paerIndiv

all: paer hsal
.PHONY: all
