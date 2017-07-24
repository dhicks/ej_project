#' ---
#' title: Company behavior and environmental justice
#' author: Daniel J. Hicks, Ph.D., <hicks.daniel.j@gmail.com>
#' output:
#'     html_document:
#'         toc: true
#' ---

#+ packages
library(tidyverse)
library(cowplot)
library(stringr)
library(broom)
library(rgdal)
library(lme4)

## Parking a memory check where it's easy to find
pryr::mem_used()

#' Load EJ data
#' --------------------
#' In this section we load the EJSCREEN dataset from EPA.  The dataset can be downloaded at <ftp://newftp.epa.gov/EJSCREEN/2016/>.  The dataset is approximately 800 MB, and contains demographic and environmental measures for 220k block groups across the US.  Descriptive names for the variables are given in the intraline comments.  
#' 
#' This section also loads a dataset that links FIPS codes to state and county names.  This dataset can be downloaded at <https://www.census.gov/geo/reference/codes/cou.html>.  
#' 
#' The plots in this section are histograms of the environmental indicators.  Most of the variables are highly skewed, and the plots indicate that log transformations make them more Gaussian.  One variable, PRE1960PCT, is a percentage; it is approximately Gaussian with a logit or log odds transformation.  
#' 
#' The JUST Capital data loaded in the next section is at the county level.  So the final step in this section is to construct county-level aggregates of the environmental data, using population-weighted averages.  Since all of the environmental indicators are log transformed, these amount to population-weighted geometic means of the un-transformed variables (where the odds is considered the un-transformed version of PRE1960PCT). 
#' 
#+ ej_load, cache = TRUE
dataf_ = read_csv('data/EJ Screen/EJSCREEN_Full_V3_USPR_TSDFupdate.csv', 
                  na = c('', 'NA', 'None'))

## Parse FIPS codes
fips_codes = read_csv('data/national_county.csv')
fips_parsed = str_match(dataf_$ID, 
                        '([0-9]{2})([0-9]{3})([0-9]{6})([0-9])') %>%
    as_tibble() %>%
    set_names(c('fips', 'state_fp', 'county_fp', 
                'tract', 'block_group')) %>%
    left_join(fips_codes)

## Variables of interest
identifiers = c('fips', 'state_fp', 'county_fp', 'tract', 
                'block_group', 'state', 'county')
ej_indicators = c('CANCER',  ## NATA Air Toxics Cancer Risk
               'RESP',       ## NATA Respiratory Hazard Index
               'DSLPM',      ## NATA Diesel Particulate Matter
               'PM25',       ## PM 2.5 Concentration Score
               'OZONE',      ## Ozone Concentration Score
               'PRE1960PCT', ## Pct. Housing Units Built Prior to 1960 (lead paint indicator)
               'PTRAF',      ## Traffic Proximity
               'PRMP',       ## RMP Proximity
               ## potential chemical accident management plan
               'PTSDF',      ## TSDF Proximity
               ## hazardous waste management facilities
               'PNPL',       ## Superfund Proximity
               'PWDIS'       ## Water Discharger Proximity
)
demographics = c('ACSTOTPOP', ## Total Population
                 'MINORPCT',  ## Pct. Minority Population
                 'LOWINCPCT', ## Pct. Low Income (<2x poverty level)
                 'LESSHSPCT', ## Less than High School Education
                 'LINGISOPCT',## Pct. Linguistically Isolated
                 'UNDER5PCT', ## Pct. Under Age 5
                 'OVER64PCT'  ## Pct. Over Age 64
)

## Normal or logged variables?  
## Normal
dataf_ %>%
    select(one_of(ej_indicators)) %>%
    gather(indicator, value) %>%
    ggplot(aes(value)) + 
    geom_histogram() +
    facet_wrap(~ indicator, scales = 'free')
## Logged
dataf_ %>%
    select(one_of(ej_indicators)) %>%
    gather(indicator, value) %>%
    ggplot(aes(log1p(value))) + 
    geom_histogram() +
    facet_wrap(~ indicator, scales = 'free')
## OZONE, PM25 work as normal, but also logged
## Everything else except PRE1960PCT works at least a little better logged
## PRE1960PCT is a percentage; it's ~normal with a logit transform
dataf_$PRE1960PCT %>%
    {log(./(1-.))} %>%
    as_tibble() %>%
    ggplot(aes(value)) + geom_histogram()

## Join EJ and FIPS datasets and select down to variables of interest
bg_dataf = dataf_ %>%
    ## Log the EJ indicators
    mutate_at(vars(one_of(ej_indicators), -PRE1960PCT), 
              funs(log10(.+1))) %>%
    ## ... except for PRE1960PCT, which gets a logit
    mutate(PRE1960PCT = log10(PRE1960PCT / (1 - PRE1960PCT) + 1)) %>%
    ## Remove infinite values in logit PRE1960PCT
    mutate(PRE1960PCT = ifelse(is.infinite(PRE1960PCT), 
                               max(PRE1960PCT[is.finite(PRE1960PCT)]), 
                               PRE1960PCT)) %>%
    bind_cols(fips_parsed) %>%
    # filter(state %in% c('DC', 'MD', 'VA'))
    select(one_of(identifiers, ej_indicators, demographics))

## Aggregate to county-level data
## NB We'll use population-weighted means, but this isn't necessarily the best option
counties_df = bg_dataf %>%
    mutate(fips = str_c(state_fp, county_fp)) %>%
    group_by(state_fp, state, county, fips) %>%
    # summarize_at(one_of(indicators), mean, na.rm = TRUE) %>% 
    summarize_at(vars(one_of(ej_indicators, demographics), 
                      -ACSTOTPOP), 
                 funs(weighted.mean(., ACSTOTPOP, na.rm = TRUE))) %>%
    ungroup()
## Need to handle population totals separately
counties_df = bg_dataf %>%
    mutate(fips = str_c(state_fp, county_fp)) %>%
    group_by(fips) %>%
    summarize(pop = sum(ACSTOTPOP)) %>%
    left_join(counties_df)
## For faster loading, discard dataf_ and cache the working datasets
rm(dataf_)


#' Load JUST Capital data
#' --------------------
#' This section loads the `COMPANY_DATA` set provided by JUST Capital.  The Excel version of this sheet was password protected and had non-descriptive variable names; a separate codebook provided slightly more informative variables names.  For convenience, a CSV was prepared manually that was not password protected and used the more-informative variable names.  However, no further details were available about these variables, including units, the way they were calculated, and their provenance.  
#' 
#' The dataset includes county-level information for 107 companies.  The company ID used in this dataset was unique, and so these data could not be linked to other JUST Capital datasets.  To account for population differences between counties, variables were normalized by the number of employees that a company has in a given county.  Log transformations were also applied.  
#' 
#+ just_load
just_df = read_csv('data/JUST Capital/COMPANY_DATA.csv') %>% 
    ## Fix deletion of leading 0 in FIPS code
    mutate(fips = str_pad(state_county_fips, 5, 
                          'left', pad = '0')) %>%
    select(fips, everything(), -state_county_fips) %>%
    ## Give company ID and number of employees more perspicuous names
    rename(company_id = id, 
           n_employees = weight1) %>%
    ## Parse currency columns
    mutate_if(is.character, 
              funs(str_replace_all(., '[\\$,]', ''))) %>%
    mutate_at(vars(-fips, -company_id, -n_employees),
              as.numeric)

just_indicators = just_df %>%
    select(-fips, -company_id, -n_employees, -living_wage) %>%
    names()

## Again, normal or logged variables? 
## Normal
just_df %>%
    select(n_employees, one_of(just_indicators)) %>%
    gather(indicator, value, -n_employees) %>%
    ggplot(aes(value/n_employees)) + 
    geom_histogram() + 
    facet_wrap(~ indicator, scales = 'free')
## Logged
just_df %>%
    select(n_employees, one_of(just_indicators)) %>%
    gather(indicator, value, -n_employees) %>%
    ggplot(aes(log1p(value/n_employees))) + 
    geom_histogram() + 
    facet_wrap(~ indicator, scales = 'free')

## Use logged
just_df = just_df %>%
    mutate_at(vars(one_of(just_indicators)), 
              funs(log10((. / n_employees)+1)))


#' Identify focal correlates
#' --------------------
#' This section identifies notable correlations between EJ and JUST Capital indicators.  Correlations are calculated between all pairs of EJ and JUST Capital indicators; in the heat map below, hotter colors (orange and red) indicate stronger absolute correlations.  The rest of the analysis focuses on the four JUST indicators that are more strongly correlated with DSLPM (which measures diesel emissions) and PTRAF (which measures traffic proximity).  
#' 
#+ id_correlates, fig.width = 10
## Join county and JUST datasets, calculate correlations, heatmap
inner_join(counties_df, just_df) %>%
    select(one_of(ej_indicators, just_indicators)) %>%
    cor(use = 'pairwise.complete.obs') %>%
    .[just_indicators, ej_indicators] %>%
    abs() %>%
    gplots::heatmap.2(scale = 'none', 
                      trace = 'none',
                      margins = c(10, 20),
                      col = gplots::rich.colors(25))
## Here we'll focus on one cluster with relatively high correlations
ej_focus = c('DSLPM',     ## Diesel particulate matter
             'PTRAF'      ## Traffic
)
just_focus = c('subsidy_sum_by_emp_total_sum',
               'subsidy_by_emp_total_sum',
               'subsidy_abs_by_emp_eitc_avg',
               'subsidy_by_emp_eitc_avg',
               'subsidy_total_eitc_sum', 
               'hourly_wage_raw_avg')

summary(counties_df[ej_focus])
summary(just_df[just_focus])


#' Plot focal correlates
#' --------------------
#' This section plots the relationship between the focal EJ and JUST variables.  The first plot shows trendlines (linear regressions) for all companies combined.  These trends show a consistent negative relationship between the EJ and JUST variable.  Since this is on a log-log scale, there is an exponential decay between increasing EJ indicator value (decreasing environmental quality) and JUST indicator value.  This suggests that, on the whole, the companies in the dataset treat their workers differently depending on the background environmental conditions (diesel emissions and vehicle traffic).  This in turn suggests a pattern of environmental discrimination; however, that depends on how the JUST variables in question should be interpreted.  
#' 
#'  The second plot shows trendlines for each individual company.  This plot shows the same general negative relationship between EJ and JUST indicators, but also shows that this relationship varies substantially across companies.  For some companies and some EJ-JUST combinations, there is a strong negative relationship; but for other companies the relationship can be basically flat.  This suggests that the degree of environmental discrimination may vary from company to company; however, again, that depends on how the JUST variables in question should be interpreted.  
#'   
#+ plot_correlates, fig.width = 16, fig.height = 10
just_counties_df = inner_join(counties_df, just_df)

## Aggregated across all companies
just_counties_df %>%
    select(company_id, fips, one_of(ej_focus, just_focus)) %>%
    gather(ej_indicator, ej_value, one_of(ej_focus)) %>%
    gather(just_indicator, just_value, one_of(just_focus)) %>%
    ggplot(aes(ej_value, just_value)) + 
    stat_binhex(aes(color = ..count..), bins = 20) +
    scale_fill_gradient(trans = 'log10') +
    geom_smooth(method = 'lm', color = 'red') +
    ylim(0, NA) +
    facet_wrap(~ ej_indicator + just_indicator, scales = 'free', 
               ncol = length(just_focus))
## Just one pair of indicators
aggregate_plot = ggplot(just_counties_df, 
                        aes(PTRAF, hourly_wage_raw_avg)) +
    geom_hex(aes(alpha = ..count..), bins = 50) +
    scale_fill_gradient(trans = 'log10', high = 'red') +
    scale_alpha_continuous(trans = 'log10') +
    geom_smooth(method = 'lm', color = 'blue') +
    ylim(0, NA)
    
## Splitting companies
just_counties_df %>%
    select(company_id, fips, one_of(ej_focus, just_focus)) %>%
    gather(ej_indicator, ej_value, one_of(ej_focus)) %>%
    gather(just_indicator, just_value, one_of(just_focus)) %>%
    ggplot(aes(ej_value, just_value)) + 
    geom_smooth(aes(group = company_id), method = 'lm') +
    facet_wrap(~ ej_indicator + just_indicator, scales = 'free', 
               ncol = length(just_focus))
## Just one pair
pair_plot = ggplot(just_counties_df, aes(PTRAF, hourly_wage_raw_avg)) +
    geom_line(stat = 'smooth', method = 'lm', se = FALSE,
                color = 'blue', alpha = .5,
                aes(group = company_id)) +
    ylim(0, NA)

## Plot the 'just one pair' plots together
# png(filename = './plot_1.png',
#     width = 8, height = 4.5, units = 'in', res = 300)
# plot_grid(aggregate_plot + #ylim(0, .6) + xlim(0, 1) +
#               xlab('Vehicle pollution') +
#               theme_cowplot(font_size = 10) +
#               theme(legend.justification=c(1,1),
#                     legend.position=c(1,1)) +
#               ggtitle('Companies combined'),
#           pair_plot + #ylim(0, .6) + xlim(0, 1) +
#               xlab('Vehicle pollution') +
#               theme_cowplot(font_size = 10) +
#               ggtitle('Companies individually'),
#           align = 'h') +
#     ggtitle('Vehicle pollution and hourly_wage_raw_avg')
# dev.off()



#' Random effects model
#' --------------------
#' The differences between companies identified in the last section can be investigated more rigorously using a random effects model  As in the second plot, this method constructs linear regressions between JUST and EJ variables for each individual company.  However, where the previous plot treated all of these regressions independently, a random effects model treats these regressions as correlated.  This helps avoid uncertainty inflation problems — we have 1 model rather than 107 independent models — and can be extended in future work, e.g., to incorporate multiple EJ and JUST indicators and control for county-level demographics.  
#' 
#' The plot shows the company-level association between vehicle pollution and the variable `hourly_wage_raw_avg`; these associations correspond to the slopes in the second plot in the last section.  The associations are expressed as the percentage association with an order-of-magnitude change in vehicle pollution.  For example, -20% means that a 10x increase in vehicle pollution corresponds to a 20% decrease in the value of `hourly_wage_raw_avg`.  The vertical bars are 95% confidence intervals on the estimates.  
#' 
#' At the most extreme, on the left end of the plot, some companies see drops of -30% or more.  
#' 
#+ rand_effects, fig.width = 10
fit = lmer(hourly_wage_raw_avg ~ 0 + (1 + PTRAF | company_id), 
     data = just_counties_df)
summary(fit)
coefs = tidy(fit, effects = 'ran_modes', conf.int = TRUE) %>%
    as_tibble()

# png(filename = './plot_2.png',
#     width = 8, height = 4, units = 'in', res = 300)
coefs %>%
    filter(term == 'PTRAF') %>%
    ggplot(aes(reorder(level, estimate), 10^estimate - 1)) + 
    geom_pointrange(aes(ymin = 10^conf.low - 1, 
                        ymax = 10^conf.high - 1), 
                    size = .25) +
    geom_hline(yintercept = 0) + 
    xlab('company ID') +
    scale_y_continuous(name = 'change per order of magnitude',
                       labels = scales::percent) +
    ggtitle('Relationship between vehicle pollution and hourly_wage_raw_avg') +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(size = 3))
# dev.off()

coefs %>%
    filter(term == 'PTRAF') %>%
    transmute(company_id = level,
              estimate = round(10^estimate - 1, digits = 3)*100) %>%
    arrange(estimate) %>%
    DT::datatable(width = 300, height = 600, rownames = FALSE)


#' Using shapefiles to plot maps
#' --------------------
#' This section contains some experimental code to plot maps.  This code was used to explore the EJSCREEN data has not yet been incorporated into this analysis.  
#' 
## NB
## census tract shapefiles are only available state by state:
## <https://www.census.gov/geo/maps-data/data/cbf/cbf_tracts.html>
## They have the format <cb_2016_ss_tract_500k.zip>
## where `ss` is the 2-digit states FIPS code
## FTP folder here:  <https://www2.census.gov/geo/tiger/GENZ2016/shp/>
# 
# ## Build mapdata for mapping
counties_shp = readOGR(dsn="data/cb_2016_us_county_500k",
                       layer="cb_2016_us_county_500k")
counties_mapdata = counties_shp %>%
    .@data %>%
    mutate(fips = as.character(GEOID)) %>%
    inner_join(fortify(counties_shp, region = 'GEOID'),
               by = c('GEOID' = 'id'))


## Mapping example
## Company 48's hourly_wage_raw_avg in California
to_map_df = just_counties_df %>%
    filter(company_id == 48, state == 'CA') %>%
    inner_join(counties_mapdata)
ggplot(to_map_df, aes(long, lat, group = group, 
                      fill = hourly_wage_raw_avg)) +
    geom_polygon()

sessionInfo()
