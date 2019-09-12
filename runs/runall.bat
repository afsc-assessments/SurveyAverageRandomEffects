set mod=rw_ln
::call run sims_wp_2013
::move sims_wp_2013.rep ln_sims_wp_2013.rep
::call run sims_pop_2013
::move sims_pop_2013.rep ln_sims_pop_2013.rep
:: --- do normal error versions
set mod=rw
call run sims_wp_2013
call run sims_pop_2013
