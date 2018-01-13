/// 786 

/******************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

#include "common.h"
#include "search.h"
#include "search_main.h"

using namespace std;

/******************************************************************************/

extern int64_t TOTAL_ATTEMPTED;
extern int64_t JACCARD_FAILED;
extern int64_t QGRAM_NORMAL_FAILED;
extern int64_t OTHER_FAILED;
extern int64_t CORE_FAILED;
extern int64_t INTERVAL_FAILED;

/******************************************************************************/

void search_main(string ref_path, string query_chr, string ref_chr, bool is_ref_complement)
{
	FastaReference fr(ref_path);

	string ref = fr.get_sequence(ref_chr);	
	auto ref_hash = make_shared<Index>(make_shared<Sequence>(ref_chr, ref, is_ref_complement));

	auto query_hash = ref_hash; 
	bool is_same_genome = (ref_chr == query_chr) && !is_ref_complement;
	if (!is_same_genome) {
		string query = fr.get_sequence(query_chr);
		query_hash = make_shared<Index>(make_shared<Sequence>(query_chr, query));
	}

	eprn("Same genome:        {}", is_same_genome);
	eprn("Reverse complement: {}", is_ref_complement);

	Tree tree;
	int total = 0, track = 0;
	for (int qi = 0; qi < query_hash->minimizers.size(); qi++) {
		auto &qm = query_hash->minimizers[qi];
		if (qm.hash.status != Hash::Status::HAS_UPPERCASE) 
			continue; // ignore N or lowercase hashes

		if (qm.loc / 10000 != track) {
			eprnn("\r   {} {:.1f}% (loci={:n} hits={:n})", 
				string(int(pct(qm.loc, query_hash->seq->seq.size()) / 2) + 1, '-'), 
				pct(qm.loc, query_hash->seq->seq.size()), qm.loc, total
			);
			track = qm.loc / 10000;
		}

		auto hits = search(qi, *ref_hash, *query_hash, tree, is_same_genome);
		for (auto &pp: hits) {
			prn("{}", pp.to_bed());
		}
		total += hits.size();
	}

	eprn("");
	eprn("Total:           = {:10n}", total);
	eprn("Fails: attempts  = {:10n}\n"
	     "       Jaccard   = {:10n}\n"
	     "       interval  = {:10n}\n"
	     "       lowercase = {:10n}\n"
	     "       q-grams   = {:10n}\n"
	     "       cores     = {:10n}",
	     TOTAL_ATTEMPTED, 
	     JACCARD_FAILED, 
	     INTERVAL_FAILED, 
	     OTHER_FAILED,
	     QGRAM_NORMAL_FAILED, 
	     CORE_FAILED
	);

	// if (query_hash != &ref_hash) {
	// 	delete query_hash;
	// }
}
