models = mnull mbatch mfull

# core recipe
sampling/hsalinarum/%/posterior_0.pkl:
	cd sampling;\
	mkdir -p hsalinarum/$*/;\
	python hsalinarum.py $(word 3, $(subst /, ,$(@D))) $(word 4, $(subst /, ,$(@D))) \
	--adapt_delta=0.95 --max_treedepth=20 \
	> hsalinarum/$*/log.out \
	2> hsalinarum/$*/log.err

# H. salinarum
conditions := standard low hi

hsal: $(foreach cond,$(conditions),$(foreach model,$(models),sampling/hsalinarum/$(cond)/$(model)/samples/posterior_0.pkl))
