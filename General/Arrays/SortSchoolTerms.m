function terms_sorted = SortSchoolTerms( terms )

 num_terms = length(terms);

seasons = cell(num_terms, 1);
years = cell(num_terms, 1);

for k = 1:num_terms
    
    words = split(terms{k}, ' ');
    
    seasons{k} = words{1};
    years{k} = words{2};
    
end

unique_years = unique(years); num_unique_years = length(unique_years);
unique_seasons = unique(seasons); num_unique_seasons = length(unique_seasons);

desired_season_order = {'Winter', 'Spring', 'Summer', 'Fall'};

terms_sorted = cell(num_terms, 1);

num_terms_sorted = 0;

for k1 = 1:num_unique_years
    
    terms_this_year = terms( strcmp(years, unique_years{k1}) ); num_terms_this_year = length(terms_this_year);
    
    if num_terms_this_year > 1
        
        seasons_this_year = seasons( strcmp(years, unique_years{k1}) );
        
        locs = GetStringSortingLocations(seasons_this_year, desired_season_order);
        
        temp_terms = terms_this_year(locs);
        
        terms_sorted(num_terms_sorted + 1:num_terms_sorted + num_terms_this_year) = temp_terms;
        num_terms_sorted = num_terms_sorted + num_terms_this_year;
        
    else
        
        terms_sorted(num_terms_sorted + 1) = terms_this_year;
        num_terms_sorted = num_terms_sorted + 1;
    end
end


end

