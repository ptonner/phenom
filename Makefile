models = mnull mbatch mfull

# core recipes
sampling/hsalinarum/combined/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p hsalinarum/$*/;\
	python hsalinarum.py $(word 3, $(subst /, ,$(@D))) $(word 4, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=20 \
	> hsalinarum/$*/log.out \
	2> hsalinarum/$*/log.err

sampling/hsalinarum/individual/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p hsalinarum/$*/;\
	python hsalinarum.py $(word 3, $(subst /, ,$(@D))) $(word 4, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=20 \
	> hsalinarum/$*/log.out \
	2> hsalinarum/$*/log.err

# H. salinarum
conditions := low hi
hbatches := $(basename $(wildcard data/hi-oxidative/*))
# hbatches := "20150517 PQ 3"  "20150630 PQ 5"  "20150704 PQ 7"  "20150717 PQ 9"     "20161107_PQ_osmo_combo" "20150607 PQ 4"  "20150702 PQ 6"  "20150715 PQ 8"  "20161010_PQ_osmo"
hbatches := 20150517\\ PQ\\ 3

hsalCombined: $(foreach cond,$(conditions),$(foreach model,$(models),sampling/hsalinarum/combined/$(cond)/$(model)/samples/posterior_0.pkl))
hsalIndiv: $(foreach cond,$(conditions),$(foreach batch,$(hbatches),sampling/hsalinarum/individual/$(cond)/$(batch)/samples/posterior_0.pkl))
hsal: hsalCombined hsalIndiv
.PHONY: hsalCombined hsalIndiv hsal
