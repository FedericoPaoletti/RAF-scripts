def ReduceToFgenerated(RPsets, final_closure):
	for reactiontuple in RPsets:
		if not reactiontuple[0].issubset(final_closure):
			RPsets.remove(reactiontuple)
	return RPsets
			
