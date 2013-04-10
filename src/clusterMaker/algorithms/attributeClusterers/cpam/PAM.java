//NB all parameters are passed by pointer
//TODO find out which variables are arrays
//TODO find out which variables are modified by function calls!

import java.io.System;


public class PAM
{

	/**
	 * Cluster by PAM.
	 * All arrays (in or out) must be preallocated to appropriate size
	 */
	public void cluster(
		// in: nrow(x); number of nodes (n)
		int nn,
		// in: ncol(x); number of attributes
		int p,
		// in: number of clusters
		int kk,
		// in: data
		double[] x,
		// out: distances, double^(nn*(nn-1)/2 + 1)
		double[] dys,
		// in:
		// 0 => compute distances from x
		// 1 => distances provided in x
		int jdyss,
		// in
		double valmd,
		// in
		int jtmd,
		// in: distance metric, 1 => euclidean; otherwise => manhattan
		int ndyst,
		// out: integer^n
		int[] nsend,
		// out: whether node is a medoid, logical^n
		int[] nrepr,
		// out: if cluster.only integer^1 else integer^n
		int[] nelem,
		// out: double^n
		double[] radus,
		// out: double^n
		double[] damer,
		// out: double^n
		double[] avsyl,
		// out: double^n
		double[] separ,
		// out: pass single value by reference
		double[] ttsyl,
		// reused variable
		// in: [cluster.only, trace.level]
		// out: values of objective function after build and after swap phases
		double[] obj,
		// inout: initial medoids (modified unless cluster.only == true)
		// medoid labels use 1-based indexing
		int[] med,
		// out: cluster labels
		int[] ncluv,
		// out: cluster information (k x 5)
		double[] clusinf,
		// out: silouette information (n x 4)
		double sylinf,
		// inout: if cluster.only integer else integer^k
		// nisol[0] != 0 indicates do.swap = true
		int[] nisol,
		// in: not used (pamonce == 1 or pamonce == 2 is not implemented)
		int pamonce
	) {

	// row dimension of cluster information matrix (ncol = 5)
	int clusinf_dim1 = kk;

	// local variables
	// if false, only do clustering and return ncluv[]
	boolean all_stats = (obj[0] == 0);
	// if true med[] contain initial medoids
	boolean med_given = (med != null);
	boolean do_swap = (nisol[0] != 0);

	int k, i, nhalf;
	int[] jhalt = new int[1];
	int trace_lev = (int) obj[1];

	// function body
	// nhalf := #{distances}+1 = sizeof(dys)
	// number of edges plus one
	nhalf = nn * (nn -1) / 2 + 1;

	if (jdyss != 1) {
		// compute distances from x

		jhalt[0] = 0;

		if (trace_lev) {
			System.out.printf("C pam(): computing %d dissimilarities: ", nhalf);
		}

		dysta(nn, p, x, dys, ndyst, jtmd, valmd, jhalt);

		if (trace_lev) {
			Sytem.out.printf("[OK]\n");
		}

		if (jhalt[0] != 0) {
			// error occured during distance computation
			// DIFF: jdyss is not set
			//jdyss = -1;
			return;
		}
	}


	double s = 0.0;
	// s := max( dys[.] ), the largest distance
	// self distance (dys[0] == 0) is not considered
	for (i = 1; i < nhalf; ++i) {
		if (s < dys[i]) s = dys[i];
	}

	// FIXME: work with med[] = (i_1, i_2, ..., i_k)
	//        instead nrepr[] = (b_1, ..., b_n)  b_i in {0,1}

	// initialize nrepr (whether node is a medoid)
	for (i = 0; i < nn; ++i) {
		nrepr[i] = 0;
	}

	// if true, med[] contain initial medoids
	if (med_given) {
		// for the moment, translate these to nrepr[] 0/1 :
		// not assuming that the med[] indices are sorted
		for (k = 0; k < kk; ++k) {
			nrepr[med[k] - 1] = 1;
		}
	}

	// Build + Swap
	// but no build if (med_given); swap only if (do_swap)
	
	bswap(kk, nn, nrepr, med_given, do_swap, trace_lev,
		radus, damer, avsyl, dys, s, obj, pamonce);
	// NB  radus, damer, and avsyl are used for scratch space here

	if (trace_lev) {
		System.out.printf("end{bsawp()}, ");
	}

	// Compute clusering & STATS if (all_stats)
	cstat(kk, nn, nsend, nrepr, all_stats,
		radus, damer, avsyl, separ, s, dys, ncluv, nelem, med, nisol);
	
	if (trace_lev) {
		System.out.printf("end{cstat()}\n");
	}

	if (all_stats) {

		// assemble cluster information into one array (matrix)
		for (k = 0; k < kk; ++k) {
			clusinf[k] = (double) nrepr[k];
			clusinf[k + clusinf_dim1] = radus[k];
			clusinf[k + (clusinf_dim1 * 2)] = avsyl[k];
			clusinf[k + (clusinf_dim1 * 3)] = damer[k];
			clusinf[k + (clusinf_dim1 * 4)] = separ[k];
		}

		if (1 < kk && kk < nn) {
			// Compute silhouette info
			dark(kk, nn, ncluv, nsend, nelem, nrepr,
				radus, damer, avsyl, ttsyl, dys, s, sylinf);
		}

	}

	/**
	 * Clustering algorithm in 2 parts: build, swap
	 *
	 * 0-based indexing
	 */
	public void bswap(
		// in: number of clusters
		int kk,
		// in: number of nodes
		int n,
		// nrepr[]: here is boolean (0/1): 1 = "is a representative object (medoid)"
		int[] nrepr,
		// in: skip build step if initial medoids are given
		boolean med_given,
		// in: skip swap unless swap is true
		boolean do_swap,
		// debugging trace level
		int trace_lev,
		// out
		// dysma[j] := D_j
		//          := d(j, medoid_nearest)
		// [KR p.102, 104]
		double[] dysma,
		// out
		// dysmb[j] := E_j
		//          := d(j, medoid_second_nearest)
		// [KR p.103]
		double[] dysmb,
		// out
		double[] beter,
		// in: distances
		double[] dys,
		// in: largest value in dys
		double s,
		// out: values of objective function before and after swap (size 2)
		double[] obj,
		// in
		int pamonce
	) {

		int i, j, ij, k, h;
		double sky;

		if (trace_lev) {
			System.out.printf("pam()'s bswap(*, s=%g, pamonce=%d): ", s, pamonce);
		}

		// ensure that s is larger than dys[i] for all i
		// (alternative is to set s to DBL_MAX, which is too large)
		s = s * 1.1 + 1.0;

		// PROPOSAL: when n is large compared to k (k == kk)
		// use a sparse representation
		// instead of boolean vector nrepr[], use ind_repr <- which(nrepr)
		
		// initialize dysma to max value
		for (i = 0; i < n; ++i) {
			dysma[i] = s;
		}

		if (med_given) {
			if (trace_lev) System.out.printf("medoids given\n");

			// dysma[j] := D(j, medoid_nearest)
			for (i = 0; i < n; ++i) {
				if (nrepr[i] == 1) {
					for (j = 0; j < n; ++j) {
						ij = ind_2(i, j);
						if (dysma[j] > dys[ij]) {
							dysma[j] = dys[ij];
						}
					}
				}
			}
		} else {
			
			// BUILD phase
			
			if (trace_lev) System.out.printf("build %d medoids:\n", kk);

			// find kk medoids
			for (k = 1; k <= kk; ++k) {
				// OPTIONAL: check user interrupt

				// compute beter[i] for all non-medoids
				// also find:
				// ammax := max_i beter[i]
				// nmax := argmax_i beter[i]

				int nmax = -1;
				double ammax, cmd;
				ammax = 0.0;

				for (i = 0; i < n; ++i) {
					if (nrepr[i] == 0) {
						// node i is not a medoid
						beter[i] = 0.0;
						for (j = 0; j < n; ++j) {
							cmd = dysma[j] - dys[ind_2(i, j)];
							if (cmd > 0.0) {
								// node i is closer to j than j's nearest medoid is to j
								beter[i] += cmd;
							}
						}
						if (ammax <= beter[i]) {
							// need "<="; "<" does not work
							// (makes a difference for ties later on...)
							ammax = beter[i];
							nmax = i;
						}
					}
				}

				// node nmax will be the new medoid
				// this non-medoid maximizes the improvement (decrease) in dysma[.]
				// (if a tie exists, the node with the largest index is chosen)
				nrepr[nmax] = 1;

				if (trace_lev >= 2) System.out.printf("    new repr. %d\n", nmax);
				
				// update dysma to reflect that node nmax is now a medoid
				for (j = 0; j < n; ++j) {
					ij = ind_2(nmax, j);
					if (dysma[j] > dys[ij]) {
						dysma[j] = dys[ij];
					}
				}

			}
			
		} // if (med_given)

		// output of above loop: nrepr[], dysma[]
		
		// DIFF: trace of nrepr and dysma not translated
		
		//dig_n = 1;


		// compute the sum of D(j, medoid_nearest)
		sky = 0.0;
		for (j = 0; j < n; ++j) {
			sky += dysma[j];
		}
		// objective function after build phase
		obj[0] = sky / n;


		// SWAP phase

		if (do_swap && (kk > 1 || med_given)) {
			
			double dzsky;
			int hbest = -1, ibest = -1, kbest = -1;
			int[] medoids = null;
			int[] clustmembership = null;
			double[] fvect = null;

			// NB: only pamonce == FALSE is translated

			// in the following, we RE-compute dysma[];
			// don't need to do so in the first iteration;
			// only need to update after swap

			while (true) {
			
				// compute dysma and dysmb
				// find distances of each node to nearest medoid (dysma) and
				// to second-nearest medoid (dysmb)
				for (j = 0; j < n; ++j) {
					// initialize to max value
					dysma[j] = s;
					dysmb[j] = s;
					for (i = 0; i < n; ++i) {
						if (nrepr[i]) {
							ij = ind_2(i, j);
							if (dysma[j] > dys[ij]) {
								dysmb[j] = dysma[j];
								dysma[j] = dys[ij];
							} else if (dysmb[j] > dys[ij]) {
								dysmb[j] = dys[ij];
							}
						}
					}
				}

				dzsky = 1.0;
				// 1 is arbitrary; only dzsky < 0 matters in the end
				// dzsky := min_{i,j} T_{i,h}

				// iterate through all non-medoids
				for (h = 0; h < n; ++h) if (nrepr[h] == 0) {
					// OPTIONAL: check user interrupt
					
					// iterate through all medoids
					for (i = 0; i < n; ++i) if (nrepr[i]) {
						double dz = 0.0;
						// dz := T_{ih} := sum_j c_{jih}  [p.104]

						// iterate through all nodes
						for (j = 0; j < n; ++j) {
							int hj = ind_2(h, j);
							ij = ind_2(i, j);
							if (dys[ij] == dysma[j]) {
								double smaller = dysmb[j] > dys[hj] ? dys[hj] : dysmb[j];
								dz += (- dysma[j] + smaller);
							} else if (dys[hj] < dysma[j]) {
								dz += (- dysma[j] + dys[hj]);
							}
						}

						if (dzsky > dz) {
							dzsky = dz;
							hbest = h;
							ibest = i;
						}
					}

				}

				final static double DBL_EPSILON = 2.2204e-16;
				if (dzsky < -16*DBL_EPSILON * fabs(sky)) {
					// essentially, dzsky < 0
					// but 'dzsky < 0' gave infinite loop:
					// swapping identical objects improved objective function (sky)

					if (trace_lev >= 2) {
						System.out.printf("    swap new %d <- %d old; decreasing diss. %7g by %g\n", hbest, ibest, sky, dzsky);
					}

					// swap non-medoid h with medoid n
					nrepr[hbest] = 1;
					nrepr[ibest] = 0;

					// update objective function
					sky += dzsky;	

				} else {
					// no improvement: break loop
					break;
				}

			}

		}

		// objective function after swap phase
		obj[1] = sky / n;

	}

	// FIXME: 1-based indexing
	/**
	 * Compute STATistics (numerical output) concerning each partition
	 */
	public void cstat(
		// in
		int kk,
		// in
		int nn,
		// out: integer^k
		// assigned medoid label of each element
		int[] nsend,
		// in
		int[] nrepr,
		// in
		boolean all_stats,
		// in
		double[] radus,
		// in
		double[] damer,
		// in
		double[] avsyl,
		// out
		double[] separ,
		// in
		double s,
		// in
		double[] dys,
		// out
		int[] ncluv,
		// out
		int[] nelem,
		// out: integer^k
		int[] med,
		// out: integer^k
		int[] nisol
	) {

		// Wall
		int j, k, ja, jk, nplac, ksmal = -1

		// make ss largest than s := max_i dys[i]
		double ss = s * 1.1 + 1.0;

		// parameter adjustments
		// DAVID: for 1-based indexing
		// FIXME convert to 0-based indexing
		--nelem;
		--ncluv;
		--separ;
		--avsyl;
		--damer;
		--radus;
		--nrepr;
		--nsend;

		// nsend[j] := i, where x[i,] is the medoid to which x[j,] belongs
		for (j = 0; j < nn; ++j) {
			if (nrepr[j] == 0) {
				double dsmal = ss;
				// iterate through medoids
				for (k = 0; k < nn; ++k) {
					if (nrepr[k] == 1) {
						int kj_ = ind_2(k, j);
						if (dsmal > dys[kj_]) {
							// assign node j to closest medoid
							dsmal = dys[kj_];
							ksmal = k;
						}
					}
				}
				nsend[j] = ksmal;
			} else {
				// node j is a medoid; assign node j to itself
				nsend[j] = j;
			}
		}

		// ncluv[j] := k, the cluster number (k = 0..(kk-1))

		jk = 0;

		// current medoid label
		nplac = nsend[0];
		
		// initialize ncluv, simutaneously assign the first cluster
		for (j = 0; j < nn; ++j) {
			ncluv[j] = 0;
			if (nsend[j] == nplac) {
				ncluv[j] = jk;
			}
		}

		// iterate from the second node
		for (ja = 1; ja < nn; ++ja) {
			// set current medoid label
			nplac = nsend[ja];
			if (ncluv[nplac] == 0) {
				// current node has not been assigned a cluster label
				++jk;
				// assign cluster label to all nodes assigned to the current medoid, starting from the second node
				for (j = 1; j < nn; ++j) {
					if (nsend[j] == nplac) {
						ncluv[j] = jk;
					}
				}
				// early break: final cluster label has already been assigned
				if (jk == kk-1) break;
			}
		}

		if (all_stats) {
			// analysis of the clustering

			int numl;
			for (k = 0; k < kk; ++k) {
				int ntt = 0, m = -1;
				double ttt = 0.0;
				radus[k] = -1.0;
				// OPTIONAL: check user interrupt here
				for (j = 0; j < nn; ++j) {
					if (ncluv[j] == k) {
						double djm;
						++ntt;
						m = nsend[j];
						nelem[ntt] = j;
						djm = dys[ind_2(j, m)];
						ttt += djm;
						if (radus[k] < djm) {
							radus[k] = djm;
						}
					}
				}
				if (ntt == 0) System.out.printf("bug in cstat(): ntt = 0\n");
				avsyl[k] = ttt / ntt;
				med[k] = m;
			}

			if (kk == 1) {
				damer[0] = s;
				nrepr[0] = nn;
				return;
			}

			// else kk > 1
			
			// numl = number of L-clusters
			numl = 0;
			for (k = 0; k < kk; ++k) {
				// identification of cluster k
				// nelem = vector of object indices
				// nel = number of objects
				int nel = 0;

				// OPTIONAL check user interrupt here

				for (j = 0; j < nn; ++j) {
					if (ncluv[j] == k) {
						nelem[nel] = j;
						++nel;
					}
				}

				nrepr[k] = nel;

				if (nel == 1) {
					int nvn = nelem[0];
					damer[k] = 0.0;
					separ[k] = ss;
					for (j = 0; j < nn; ++j) {
						if (j != nvn) {
							int mevj = ind_2(nvn, j);
							if (separ[k] > dys[mevj]) {
								separ[k] = dys[mevj];
							}
						}
					}

					// Is cluster k...
					// 1) an L-cluster, or
					// 2) an L*-cluster?
					if (separ[k] == 0.0) {
						++numl;
					}
				} else {
					// nel != 1
					double dam = -1.0, sep = ss;
					boolean kand = true;
					for (ja = 0; ja < nel; ++ja) {
						int jb, nvna = nelem[ja];
						double aja = -1.0, ajb = ss;
						for (jb = 0; jb < nn; ++jb) {
							int jndz = ind_2(nvna, jb);
							if (ncluv[jb] == k) {
								if (aja < dys[jndz]) {
									aja = dys[jndz];
								}
							} else {
								if (ajb > dys[jndz]) {
									ajb = dys[jndz];
								}
							}
						}
						if (kand && aja >= ajb) {
							kand = false;
						}
						if (dam < aja) {
							dam = aja;
						}
						if (sep > ajb) {
							sep = ajb;
						}
					}
					separ[k] = sep;
					damer[k] = dam;
					if (kand) {
						++numl;
						if (dam >= sep) {
							nisol[k] = 1;
						} else {
							nisol[k] = 2;
							continue;
						}
					}
				}
				// nel = 1 or (!kand)
				nisol[k] = 0;

			} // for (k)

		} // all_stats

	}

	// FIXME: array slices, 1-based indexing
	/**
	 * Compute silhouette information
	 */
	public void dark(
		int kk,
		int nn,
		int[] ncluv,
		int nsend,
		int nelem,
		int negbr,
		double syl,
		double srank,
		double[] avsyl,
		double[] ttsyl,
		double[] dys,
		double s,
		double[] sylinf
	) {

		int k, nsylr;

		// pointers to sylinf[] columns -- sylinf[nn, 4]
		// sylinf_1 sylinf_2, sylinf_3, sylinf_4;
		
		// parameter adjustments
		// David FIXME: These are meant to point to one element ahead of the first element of the array....
		//              Cannot do this in Java!
		//              These adjustments were for 1-based indexing
		--avsyl;
		--ncluv;

		nsylr = 0;
		ttysl = 0.0;

		for (k = 1; k <= kk; ++k) {

			// nelem[0:(ntt-1)] := indices (1-based) of obs in cluster k
			int j, l, ntt = 0;
			for (j = 1; j <= nn; ++j) {
				if (ncluv[j] == k) {
					nelem[ntt] = j;
					++ntt;
				}
			}

			// (j+1)-th obs. in cluster k
			for (j = 0; j < ntt; ++j) {
				int k_, nj = nelem[j];
				// WHY 1.1?
				double dysb = s * 1.1 + 1.;
				negbr[j] = 1;

				// for all clusters k_ != k
				for (k_ = 1; k_ <= kk; ++k_) if (k_ != k) {
					double db = 0.0;
					int nbb = 0;
					for (l = 1; l <= nn; ++l) if (ncluv[l] == k_) {
						++nbb;
						if (l != nj) {
							// NOTE
							db += dys[ind_2(nj, l)];
						}
					}
					// now db(k_) := mean( d[j, l]; l in c_{k_} )
					db /= nbb;
					if (dysb > db) {
						dysb = db;
						negbr[j] = k_;
					}
					// negbr[j] := arg max_{k_} db(k_)
				}

				if (ntt > 1) {
					double dysa = 0.0;
					for (l = 0; l < ntt; ++l) {
						int nl = nelem[l];
						if (nj != nl) {
							dysa += dys[ind_2(nj, nl)];
						}
					}
					dysa /= ntt - 1;
					if (dysa > 0.0) {
						if (dysb > 0.0) {
							if (dysb > dysa) {
								syl[j] = 1.0 - dysa / dysb;
							} else if (dysb < dysa) {
								syl[j] = dysb / dysa - 1.0;
							} else {
								// dysb == dysa
								syl[j] = 0.0;
							}
							if (syl[j] < -1.0) {
								syl[j] = -1.0;
							} else if (syl[j] > 1.0) {
								syl[j] = 1.0;
							}
						} else {
							syl[j] = -1.0;
						}
					} else if (dysb > 0.0) {
						// dysa == 0
						syl[j] = 1.0;
					} else {
						syl[j] = 0.0;
					}
				} else {
					// ntt == 1
					syl[j] = 0.0;
				}

			} // for (j)

			avsyl[k] = 0.0;

			if (ntt == 0) {
				// this can happen when medoids are user-specified
				// next k
				continue
			}

			for (j = 0; j < ntt; ++j) {
				// Wall
				int lang = -1;
				double symax = -2.0;
				for (l = 0; l < ntt; ++l) {
					if (symax < syl[l]) {
						symax = syl[l];
						lang = l;
					}
				}
				nsend[j] = lang;
				srank[j] = symax;  // = syl[lang]
				avsyl[k] += srank[j];
				syl[lang] = -3.0;
			}

			ttsyl[0] += avsyl[k];
			avsyl[k] /= ntt;
			if (ntt == 1) {
				// FIXME convert (x,y) index to linear index
				sylinf[nsylr, 0] = (double) k;
				sylinf[nsylr, 1] = (double) negbr[0];
				sylinf[nsylr, 2] = 0.0;
				sylinf[nsylr, 3] = (double) nelem[0];
				++nsylr;
			} else {
				for (j = 0; j < ntt; ++j) {
					int lplac = nsend[j];
					sylinf[nsylr, 0] = (double) k;
					sylinf[nsylr, 1] = (double) negbr[lplac];
					sylinf[nsylr, 2] = srank[j];
					sylinf[nsylr, 3] = (double) nelem[lplac];
					++nsylr;
				}
			}

		} // for (k)

		ttsyl[0] /= nn;

	}

	// DONE
	// NB converted from 1-based index to 0-based index
	/**
	 * Compute distances array.
	 *
	 * 0-based indexing
	 */
	public void dysta(
		// in: number of nodes
		int nn,
		// in: number of attributes
		int p,
		// in: data
		double[] x,
		// out: distances, preallocated with size double^{m},
		// where m = number of edges + 1 = nn*(nn-1)/2 + 1
		// one extra space is reserved for distance to self (first element)
		double[] dys,
		// in:  1 => eucliden; otherwise => manhattan
		int ndyst,
		// in: integer^p
		// -1 => column/attribute contains an NA
		//  1 => column/attribute contains no NA
		int[] jtmd,
		// value for missing data (possibly a different value for each attribute)
		double[] valmd,
		// out: error flag (non-zero value denotes error)
		int[] jhalt
	) {

		int nlk, j, l, k, lsubt, npres;
		double pp, clk;

		// start at index 0
		nlk = 0;

		// first element denotes distance to self:
		// d[i,i] == dys[0] == 0
		dys[0] = 0.0;
		
		// iterate from the second element
		for (l = 1; l < nn; ++l) {
			// number of other nodes
			lsubt = l-1;
			// iterate from the first element
			for (k = 0; k < lsubt; ++k) {
				// distance value to be calculated
				clk = 0.0;
				// current position in distance array
				nlk += 1;
				// number of non-NA data points
				npres = 0;
				// iterate through attributes
				for (j = 0; j < p; ++j) {
					// current attribute contains at least one NA value
					// check if current data point is NA
					if (jtmd[j] < 0) {
						// exact floating point equality
						// admissible because binary representations are identical
						if (x[l, j] == valmd[j]) break;
						if (x[k, j] == valmd[j]) break;
						// data point is NA => ignore data point
					}
					npres += 1;
					double d = x[l,j] - x[k,j];
					if (ndyst == 1) {
						clk += d * d;
					} else {
						clk += Math.abs(d);
					}
				}
				if (npres == 0) {
					// all data points are NA for at least one distance calculation
					// set the error flag
					jhalt = 1;
					// arbitrarily set distance to -1
					dys[nlk] = -1.0;
				} else {
					// scale the distance to the number of non-NA data points
					// double cast is not needed,
					//   since evaluation starts from a double value (clk)
					if (ndyst == 1) {
						dys[nlk] = Math.sqrt(clk * p / npres);
					} else {
						dys[nlk] = clk * p / npres;
					}
				}
			}
		}

	}

	/**
	 * Convert 2-D coordinate into linear index of lower diagonal matrix (starting at 1,0 excluding diagonal).
	 * Linear index of dys[.] where d(i,j) is stored:
	 * d(l,j) == dys[ind_2(i,j)]
	 *
	 * 0-based indexing
	 */
	private int ind_2(int i, int j) {
		// max_m is the largest integer m s.t. (m-1)*m < Integer.MAX_VALUE = 2^31 - 1
		static final max_m = 46341;

		int result = 0;
		if (i != j) {
			// m := larger index
			int m = (i > j) ? i : j;
			// n := smaller index
			int n = (i > j) ? j : i;

			result = (m <= max_m)
				? (m-1)*m/2 + n + 1
				: (int) ((double)(m-1)*m/2 + n + 1);
		}

		return result;
	}

}

