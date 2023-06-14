COMP = c++
all: generate_matching_instances generate_compression_matchings generate_pairedmatchings join_pairedmatchings remove_equivalent_matchedpairs remove_equivalent_exhaust
generate_matching_instances: generate_matching_instances.cpp
	$(COMP) generate_matching_instances.cpp -O3 -std=c++0x -lfftw3 -o generate_matching_instances
generate_compression_matchings: generate_compression_matchings.cpp
	$(COMP) generate_compression_matchings.cpp -O3 -std=c++0x -o generate_compression_matchings
generate_pairedmatchings: generate_pairedmatchings.cpp decomps.h
	$(COMP) generate_pairedmatchings.cpp -O3 -std=c++0x -lfftw3 -o generate_pairedmatchings
join_pairedmatchings: join_pairedmatchings.cpp decomps.h
	$(COMP) join_pairedmatchings.cpp -O3 -std=c++0x -o join_pairedmatchings
remove_equivalent_matchedpairs: remove_equivalent_matchedpairs.cpp decomps.h coprimelist.h
	$(COMP) remove_equivalent_matchedpairs.cpp -O3 -std=c++0x -o remove_equivalent_matchedpairs
remove_equivalent_exhaust: remove_equivalent_exhaust.cpp decomps.h coprimelist.h
	$(COMP) remove_equivalent_exhaust.cpp -O3 -std=c++0x -o remove_equivalent_exhaust
clean:
	rm -f generate_matching_instances generate_compression_matchings generate_pairedmatchings join_pairedmatchings remove_equivalent_matchedpairs remove_equivalent_exhaust
