In this folder, we are going to store different kind of CSV files 
in order to test our code.

At the moment, our first approach is to write a code able to identify
the next types ID variants(GRCh38) :


TYPES                    e.g.

Transcript (:c.)        NM_016247.3:c.3183C>G
RefSeq Gene (:g.)       NG_028284.1:g.92745C>G
Protein (:p.)           NP_057331.2:p.(Cys1061Trp)
Protein (:p.)	        NP_057331.2:p.(C1061W)	
	


To test our code, The next 4 variants have been selected

Transcript                           RefSeq Gene (:g.)               Protein (:p.)                        Protein (:p.) 

NM_016247.3:c.3183C>G                NG_028284.1:g.92745C>G           NP_057331.2:p.(Cys1061Trp)           NP_057331.2:p.(C1061W)
NM_007294.3:c.5066T>C                NG_005905.2:g.150368T>C          NP_009225.1:p.(Met1689Thr)           NP_009225.1:p.(M1689T)
NM_000138.4:c.356G>A                 NG_008805.2:g.50564G>A           NP_000129.3:p.(Cys119Tyr)            NP_000129.3:p.(C119Y) 
NM_000257.3:c.1988G>A                N.a                              NP_000248.2:p.(Arg663His)            NP_000248.2:p.(R663H)


A 5th CSV (called ALL) has been created that store all of these variants types. 
This has been created to test the code when in a single CSV more than one
variant type is stored.
