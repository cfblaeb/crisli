def gc_filter(grnas: list, argn: int, lower: float=0.2, upper: float=0.8):
	"""
	Filter based on %GC
	:param grnas: list of grna tuples (gid, seq, strand, pos, offscore)
	:param argn:  desired number of grna
	:param lower: lowest allowed gc
	:param upper: highest allowed gc
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


def offtarget_score_filter(grnas: list, argn: int, threshold: int):
	"""
	Filter based on offtarget score
	:param grnas: list of grna tuples (gid, seq, strand, pos, offscore)
	:param argn:  desired number of grna
	:param threshold: a threshold over which targets will be filtered
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


def five_mark_filter(grnas: list, argn: int, region_start_pos: int, max_distance: int):
	"""
	Filter based on position. Position should, from 5' end (strand specific), be no more than max_distance downstream
	:param grnas: list of grna tuples (gid, seq, strand, pos, offscore)
	:param argn:  desired number of grna
	:param region_start_pos: region start bp
	:param max_distance: distance from region start (strand specific)
	:return: filtered_grnas_ids (grnas that failed this test)
	"""
	filtered_grna_ids = []


	return filtered_grna_ids