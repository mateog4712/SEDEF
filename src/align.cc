/// 786

/******************************************************************************/

#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <glob.h>

#include "align.h"
#include "common.h"
#include "fasta.h"
#include "extern/ksw2.h"

using namespace std;

/******************************************************************************/

Alignment Alignment::from_cigar(const string &a, const string &b, const string &cigar_str)
{
	auto aln = Alignment{ "A", 0, (int)a.size(), "B", 0, (int)b.size(), a, b, "", "", "", {}, {} };
	for (int ci = 0, num = 0; ci < cigar_str.size(); ci++) {
		if (isdigit(cigar_str[ci])) {
			num = 10 * num + (cigar_str[ci] - '0');
		} else if (cigar_str[ci] == ';') {
			continue;
		} else {
			aln.cigar.push_back({cigar_str[ci], num});
			num = 0;
		}
	}
	aln.populate_nice_alignment();
	return aln;
}

string Alignment::cigar_string()
{
	string res;
	for (auto &p: cigar) {
		res += fmt::format("{}{}", p.second, p.first);
	}
	return res;
}

Alignment Alignment::trim() 
{
	auto result = *this;
	while (result.cigar.size()) {
		if (result.cigar[0].first == 'D') {
			result.a = result.a.substr(result.cigar[0].second);
			result.start_a += result.cigar[0].second;
			result.cigar.pop_front();
		} else if (result.cigar[0].first == 'I') {
			result.b = result.b.substr(result.cigar[0].second);
			result.start_b += result.cigar[0].second;
			result.cigar.pop_front();
		} else if (result.cigar.back().first == 'D') {
			result.end_a -= result.cigar.back().second;
			result.a = result.a.substr(0, result.a.size() - result.cigar.back().second);
			result.cigar.pop_back();
		} else if (result.cigar.back().first == 'I') {
			result.end_b -= result.cigar.back().second;
			result.b = result.b.substr(0, result.b.size() - result.cigar.back().second);
			result.cigar.pop_back();
		} else {
			break;
		}
	}

	result.populate_nice_alignment();
	return result;
}

void Alignment::populate_nice_alignment()
{
	align_a = "";
	align_b = "";
	alignment = "";
	int ia = 0, ib = 0;
	for (auto &c: cigar) {
		for (int i = 0; i < c.second; i++) {
			assert(c.first != 'M' || ia < a.size());
			assert(c.first != 'M' || ib < b.size());
			if (c.first == 'M' && toupper(a[ia]) == toupper(b[ib])) {
				alignment += "|";
			} else {
				alignment += "*";
			}
			if (c.first != 'D') align_b += b[ib++];
			else                align_b += "-";
			if (c.first != 'I') align_a += a[ia++];
			else                align_a += "-";
		}
	}
}

Alignment::AlignmentError Alignment::calculate_error()
{
	if (!alignment.size()) populate_nice_alignment();

	int gaps = 0, gap_bases = 0, mismatches = 0, matches = 0;
	for (auto &c: cigar) {
		if (c.first != 'M') gaps++, gap_bases += c.second;
	}
	for (int i = 0; i < alignment.size(); i++) {
		if (align_a[i] != '-' && align_b[i] != '-') {
			if (toupper(align_a[i]) == toupper(align_b[i])) {
				matches++;
			} else {
				mismatches++;
			}
		}
	}
	return AlignmentError{gaps, gap_bases, mismatches, matches};
}

string Alignment::print(int width)
{
	if (!alignment.size()) {
		populate_nice_alignment();
	}

	string res;
	int qa = start_a, sa = 0;
	int qb = start_b, sb = 0;

	auto err = calculate_error();
	res += fmt::format(
		"       A: {:>9}..{:<9} (len {:7})    Gaps:       {:5} = {:.0f}% ({})\n"
		"       B: {:>9}..{:<9} (len {:7})    Mismatches: {:5} = {:.0f}%\n"
		"   CIGAR: {}\n",
		start_a, end_a, end_a - start_a, 
		err.gap_bases, err.gap_error(), err.gaps,
		start_b, end_b, end_b - start_b,
		err.mismatches, err.mis_error(),
		cigar_string()
	);
	for (int i = 0; i < alignment.size(); i += width) {
		res += fmt::format(
			"   {:10}: {} {}\n   {:10}  {} {}\n   {:10}: {} {}\n", 
			qa, align_a.substr(i, width),   sa,  
			"", alignment.substr(i, width), i+align_a.substr(i, width).size(),
			qb, align_b.substr(i, width), sb);
		for (auto c: align_a.substr(i, width)) if (c != '-') qa++, sa++;
		for (auto c: align_b.substr(i, width)) if (c != '-') qb++, sb++;
	}
	return res;
}

string Alignment::print_only_alignment(int width)
{
	if (!alignment.size()) {
		populate_nice_alignment();
	}

	if (width == -1) { 
		width = alignment.size();
	}
	string res;
	for (int i = 0; i < alignment.size(); i += width) {
		res += fmt::format(
			"{}\n{}\n{}\n\n", 
			align_a.substr(i, width),  
			alignment.substr(i, width),
			align_b.substr(i, width)
		);
	}
	return res;
}

void Alignment::cigar_from_alignment() 
{
	cigar.clear();
	int sz = 0;
	char op = 0, top;
	for (int i = 0; i < alignment.size(); i++) {
		if (align_a[i] == '-') {
			top = 'I';
		} else if (align_b[i] == '-') {
			top = 'D';
		} else {
			top = 'M';
		}

		if (op != top) {
			if (op) cigar.push_back(make_pair(op, sz));
			op = top, sz = 0;
		}
		sz++;
	}
	cigar.push_back(make_pair(op, sz));
}

vector<Alignment> Alignment::max_sum(int min_span)
{
	assert(alignment.size());

	// start/end in the alignment
	int best_score = -1, best_start = 0, best_end = 0;
	int score_so_far = 0, start_so_far = 0, end_so_far = 0;

	vector<int> lookup_a(alignment.size() + 1, 0); int ia = 0;
	vector<int> lookup_b(alignment.size() + 1, 0); int ib = 0;

	vector<pair<int, int>> hits; // start, end
	for (int i = 0; i < alignment.size(); i++) {
		int score = (alignment[i] == '|' ? 1 : -1);

		if (score > score_so_far) {
			if (best_end - best_start >= min_span 
				&& (!hits.size() || hits.back() != make_pair(best_start, best_end))) 
			{
				hits.push_back({best_start, best_end});
			}
			score_so_far = score;
			start_so_far = i;
			end_so_far = i + 1;
			best_score = score_so_far;
			best_start = start_so_far, best_end = end_so_far;
		} else {
			score_so_far += score;
			end_so_far++;
		}

		if (score_so_far > best_score) {
			best_score = score_so_far;
			best_start = start_so_far, best_end = end_so_far; 
		}

		lookup_a[i] = ia;
		lookup_b[i] = ib;
		if (align_a[i] != '-') ia++; 
		if (align_b[i] != '-') ib++;
	}
	lookup_a[alignment.size()] = ia;
	lookup_b[alignment.size()] = ib;

	if (best_end - best_start >= min_span 
		&& (!hits.size() || hits.back() != make_pair(best_start, best_end))) 
	{
		hits.push_back({best_start, best_end});
	}

	vector<Alignment> results;
	for (auto &h: hits) {
		// translate: a --> b
		results.push_back({
			chr_a, start_a + lookup_a[h.first], start_a + lookup_a[h.second],
			chr_b, start_b + lookup_b[h.first], start_b + lookup_b[h.second],
			a.substr(lookup_a[h.first], lookup_a[h.second] - lookup_a[h.first]),
			b.substr(lookup_b[h.first], lookup_b[h.second] - lookup_b[h.first]),
			align_a.substr(h.first, h.second - h.first),
			align_b.substr(h.first, h.second - h.first),
			alignment.substr(h.first, h.second - h.first),
			{}
		});
		results.back().cigar_from_alignment();
	}

	return results;
}

void Alignment::prepend_cigar(const deque<pair<char, int>> &app)
{
	assert(app.size());
	if (cigar.size() && cigar.front().first == app.back().first) {
		cigar.front().second += app.back().second;
		cigar.insert(cigar.begin(), app.begin(), app.end() - 1);
	} else {
		cigar.insert(cigar.begin(), app.begin(), app.end());
	}
}

void Alignment::append_cigar(const deque<pair<char, int>> &app)
{
	assert(app.size());
	if (cigar.size() && cigar.back().first == app.front().first) {
		cigar.back().second += app.front().second;
		cigar.insert(cigar.end(), next(app.begin()), app.end());
	} else {
		cigar.insert(cigar.end(), app.begin(), app.end());
	}
}

Alignment Alignment::from_anchors(const string &qstr, const string &rstr,
	const list<pair<int, int>> &query_kmers, const list<pair<int, int>> &ref_kmers,
	const int side) 
{
	auto qk_p = query_kmers.begin();
	auto rk_p = ref_kmers.begin();
	Alignment aln {
		"A", qk_p->first, qk_p->second, 
		"B", rk_p->first, rk_p->second,
		qstr.substr(qk_p->first, qk_p->second - qk_p->first),
		rstr.substr(rk_p->first, rk_p->second - rk_p->first),
		"", "", "",
		{{'M', qk_p->second - qk_p->first}}, {}
	};
	assert(query_kmers.size() == ref_kmers.size());

	for (auto qk = next(qk_p), rk = next(rk_p); qk != query_kmers.end(); qk++, rk++) {
		// eprnn("--- {:6}..{:6} vs {:6}..{:6}", qk_p->first, qk_p->second,
		// 	rk_p->first, rk_p->second);
		assert(qk->first >= qk_p->second);
		assert(rk->first >= rk_p->second);
		aln.end_a = qk->second;
		aln.end_b = rk->second;
		aln.a += qstr.substr(qk_p->second, qk->second - qk_p->second);
		aln.b += rstr.substr(rk_p->second, rk->second - rk_p->second);
		if (qk->first - qk_p->second && rk->first - rk_p->second) {
			auto gap = align(
				qstr.substr(qk_p->second, qk->first - qk_p->second), 
				rstr.substr(rk_p->second, rk->first - rk_p->second), 
				5, -4, 40, 1
			);
			// eprnn(" {}", gap.cigar_string());
			aln.append_cigar(gap.cigar);
		} else if (qk->first - qk_p->second) {
			aln.append_cigar({{'D', qk->first - qk_p->second}});	
		} else if (rk->first - rk_p->second) {
			aln.append_cigar({{'I', rk->first - rk_p->second}});	
		} 
		// eprn("");
		assert(qk->second - qk->first == rk->second - rk->first);
		aln.append_cigar({{'M', qk->second - qk->first}});
		qk_p = qk, rk_p = rk;
	}
	// Add end


	int qlo = query_kmers.front().first, 
		 qhi = query_kmers.back().second;
	int rlo = ref_kmers.front().first,   
		 rhi = ref_kmers.back().second;
	if (side) {
		int qlo_n = max(0, qlo - side);
		int rlo_n = max(0, rlo - side);
		if (qlo - qlo_n && rlo - rlo_n) {
			auto gap = align(
				qstr.substr(qlo_n, qlo - qlo_n), 
				rstr.substr(rlo_n, rlo - rlo_n), 
				5, -4, 40, 1
			);
			//eprn("lside={}, {}->{} / {}->{}, cigar={}", side, qlo_n, qlo, rlo_n, rlo, gap.cigar_string());
			int i = gap.cigar.size() - 1;
			int sq = 0, sr = 0;
			if (gap.cigar[i].first != 'M' && gap.cigar[i].second <= 10) {
				if (gap.cigar[i].first == 'D') {
					sq += gap.cigar[i].second;
				} else {
					sr += gap.cigar[i].second;
				}
				i--;
			}
			if (i >= 0 && gap.cigar[i].first == 'M') {
				//eprn("added!");
				sq += gap.cigar[i].second;
				sr += gap.cigar[i].second;

				aln.prepend_cigar(deque<pair<char, int>>(gap.cigar.begin() + i, gap.cigar.end()));
				aln.a = qstr.substr(qlo - sq, sq) + aln.a;
				aln.b = rstr.substr(rlo - sr, sr) + aln.b;
				aln.start_a = qlo = qlo - sq; 
				aln.start_b = rlo = rlo - sr;	
			}
		}

		int qhi_n = min(qhi + side, (int)qstr.size());
		int rhi_n = min(rhi + side, (int)rstr.size());
		if (qhi_n - qhi && rhi_n - rhi) {
			auto gap = align(
				qstr.substr(qhi, qhi_n - qhi), 
				rstr.substr(rhi, rhi_n - rhi), 
				5, -4, 40, 1
			);
			int i = 0;
			int sq = 0, sr = 0;
			//eprn("rside={}, {}->{} / {}->{}, cigar={}", side, qhi, qhi_n, rhi, rhi_n, gap.cigar_string());
			if (gap.cigar[i].first != 'M' && gap.cigar[i].second <= 10) {
				if (gap.cigar[i].first == 'D') {
					sq += gap.cigar[i].second;
				} else {
					sr += gap.cigar[i].second;
				}
				i++;
			}
			if (i < gap.cigar.size() && gap.cigar[i].first == 'M') {
				sq += gap.cigar[i].second;
				sr += gap.cigar[i].second;

				aln.append_cigar(deque<pair<char, int>>(gap.cigar.begin(), gap.cigar.begin() + i + 1));
				aln.a += qstr.substr(qhi, sq);
				aln.b += rstr.substr(rhi, sr);
				aln.end_a = qhi = qhi + sq;
				aln.end_b = rhi = rhi + sr;	
			}
		}
	}

	assert(qlo >= 0);
	assert(rlo >= 0);
	assert(qhi <= qstr.size());
	assert(rhi <= rstr.size());
	assert(aln.a == qstr.substr(qlo, qhi - qlo));
	assert(aln.b == rstr.substr(rlo, rhi - rlo));

	aln = aln.trim();
	aln.error = aln.calculate_error();
	return aln;
}

/******************************************************************************/

Alignment align_helper(const string &tseq, const string &qseq, int sc_mch, int sc_mis, int gapo, int gape, int bandwidth)
{
	const int STEP = 50 * 1000; // Max. alignment size (if larger, split into pieces)

	int8_t a = (int8_t)sc_mch, b = sc_mis < 0 ? (int8_t)sc_mis : (int8_t)(-sc_mis); // a>0 and b<0
	int8_t mat[25] = { 
		a, b, b, b, 0, 
		b, a, b, b, 0, 
		b, b, a, b, 0, 
		b, b, b, a, 0, 
		0, 0, 0, 0, 0 
	};
	deque<pair<char, int>> cigar;
	for (int SP = 0; SP < min(tseq.size(), qseq.size()); SP += STEP) {
		ksw_extz_t ez;
		ksw_extz2_sse(0, 
			min(STEP, (int)(qseq.size() - SP)), (const uint8_t*)(qseq.c_str() + SP), 
			min(STEP, (int)(tseq.size() - SP)), (const uint8_t*)(tseq.c_str() + SP), 
			5, mat, // M; MxM matrix
			gapo, gape, 
			bandwidth, -1, // band width; off-diagonal drop-off to stop extension (-1 to disable)
			0, &ez);
		for (int i = 0; i < ez.n_cigar; i++) {
			cigar.push_back({"MID"[ez.cigar[i] & 0xf], ez.cigar[i] >> 4});
		}
		free(ez.cigar);
	}
	return Alignment{ 
		"A", 0, (int)qseq.size(), 
		"B", 0, (int)tseq.size(), 
		qseq, tseq, "", "", "", cigar 
	};
}

Alignment align(const string &fa, const string &fb, int match, int mismatch, int gap_open, int gap_extend, int bandwidth)
{
	string xa = fa, xb = fb;

	transform(xa.begin(), xa.end(), xa.begin(), align_dna);
	transform(xb.begin(), xb.end(), xb.begin(), align_dna);

	auto a = align_helper(fa, fb, match, mismatch, gap_open, gap_extend, bandwidth);
	a.a = fa, a.b = fb;
	a.populate_nice_alignment();
	a.error = a.calculate_error();
	return a;
}
