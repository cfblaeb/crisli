def gc_filter(grnas: list, argn: int, lower: float=0.2, upper: float=0.8, hard_filter: bool = False):
	"""
	Filter based on %GC
	:param grnas: list of grna tuples (gid, seq, strand, pos, offscore)
	:param argn:  desired number of grna
	:param lower: lowest allowed gc
	:param upper: highest allowed gc
	:param hard_filter: If true, remove all grnas that fit criteria. If false, remove only until theres argn grna left
	:return: filtered_grnas_ids (grnas that failed this test)
	"""
	filtered_grna_ids = []  # will be filled with grna to remove from db
	# calc gc
	pgcs = [(x[1].count("C") + x[1].count("G")) / len(x[1]) for x in grnas]
	# change %gc to how far away you are from edges. E.g. 15% gc = 5% away from 20. and 85% gc is 5% away from 80...
	pgcs_edge = []
	for gc in pgcs:
		if gc < 0.2:
			pgcs_edge.append(lower - gc)
		elif gc > 0.8:
			pgcs_edge.append(gc - upper)
		else:
			pgcs_edge.append(0)
	grnas = [[*x, pgcs_edge[i]] for i, x in enumerate(grnas)]  # add percent gc to grnas
	grnas = sorted(grnas, key=lambda x: x[5], reverse=True)  # sort by how bad the gc is
	for gid, seq, strand, pos, offscore, gc_edge in grnas:
		if gc_edge != 0:
			filtered_grna_ids.append(gid)
			if len(grnas) - len(filtered_grna_ids) == argn:  # if the remaining amount of grnas is now n then stop
				break
		else:  # we have reached the end of bad gcs
			break

	return filtered_grna_ids


def offtarget_score_filter(grnas: list, argn: int, threshold: int, hard_filter: bool = False):
	"""
	Filter based on offtarget score
	:param grnas: list of grna tuples (gid, seq, strand, pos, offscore)
	:param argn:  desired number of grna
	:param threshold: a threshold over which targets will be filtered
	:param hard_filter: If true, remove all grnas that fit criteria. If false, remove only until theres argn grna left
	:return: filtered_grnas_ids (grnas that failed this test)
	"""
	filtered_grna_ids = []
	grnas = sorted(grnas, key=lambda x: x[4], reverse=True)  # sort by offscores
	for gid, seq, strand, pos, offscore in grnas:
		if offscore > threshold:
			filtered_grna_ids.append(gid)
			if len(grnas) - len(filtered_grna_ids) == argn:  # if the remaining amount of grnas is now n then stop
				break
		else:  # we have reached the end of bad gcs
			break

	return filtered_grna_ids


def five_mark_filter(grnas: list, argn: int, region_start_pos: int, region_end_pos: int, max_distance_bp: int, max_distance_pct: float, hard_filter: bool = False):
	"""
	Filter based on position. Position should, from 5' end (strand specific), be no more than max_distance downstream
	Max distance should be defined as both bp and pct. Which ever comes first will be used.
	:param grnas:
	:param argn:
	:param region_start_pos:
	:param region_strand: string of either '1' or '-1'
	:param max_distance_bp:
	:param max_distance_pct: 0-1
	:param hard_filter: If true, remove all grnas that fit criteria. If false, remove only until theres argn grna left
	:return:
	"""

	filtered_grna_ids = []
	rl = abs(region_end_pos-region_start_pos)
	if hard_filter:
		for gid, seq, strand, pos, offscore in grnas:
			dist_from_start = abs(pos - region_start_pos)
			if dist_from_start > max_distance_bp or dist_from_start/rl > max_distance_pct:
				filtered_grna_ids.append(gid)

	return filtered_grna_ids
