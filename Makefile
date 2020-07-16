# Data setup

data/hsal:
	cd data/raw/hsal; python parse.py

data/paer:
	cd data/raw/paer; python parse.py

data: data/paer data/hsal
.PHONY: data data/paer data/hsal

# Posterior inference


models = mnull mbatch mfull

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

# H. salinarum
conditions := low hi
hbatches := $(foreach batch,$(wildcard data/hi-oxidative/*),$(lastword $(subst /, ,$(batch))))

hsalCombined: $(foreach cond,$(conditions),$(foreach model,$(models),sampling/hsalinarum/combined/$(cond)/$(model)/samples/posterior_0.pkl))
hsalIndiv: $(foreach cond,$(conditions),$(foreach batch,$(hbatches),sampling/hsalinarum/individual/$(cond)/$(batch)/samples/posterior_0.pkl))
hsal: hsalCombined hsalIndiv
.PHONY: hsalCombined hsalIndiv hsal
