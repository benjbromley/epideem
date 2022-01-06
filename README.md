# epideem
Code to display and model state-level COVID-19 data

Runs under python3.6

* Fits SIR models to "wild" (alpha) and delta-variant COVID-19 case counts as reported 
   in the nytimes/covid-19-data repo. To mitigate reporting issues a 7-day avg is applied.
* Adds an omicron variant model derived from delta-variant fits boosted in both infection 
  and recovery rates by a fact of 1.2-to-1.3 (depending on version) and 1.1, respectively.
  The starting point of the omicron rise is set by 7-day avg case counts on 2022-01-01
* Generates figure "tmp.png" 

Usage: epi2022.py [State]
       Command line argument State is optional name of U.S. state (one argument, capitalized). 
       Examples: Utah, "New York"
       Default: Utah
       
 Notes:
 
 Future upgrade will perform actual SIR fit to data around 2022-01-01 as the omicron variant
 dominates over delta.
