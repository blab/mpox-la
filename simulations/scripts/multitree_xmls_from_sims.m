% getSeqFromSims

% defines the start of the simulation
sample_cutoff = {'2023-01-14'};
end_date = '2022-05-02';

% define the reporting delay in days
reporting_delay = 0;

scheme = {'half', 'all'};

for v = 1: length(scheme)
    for rep = 2:8
        for sc = 1 : length(sample_cutoff)
            rate_shifts = [7/366:7/366:(datenum(sample_cutoff(sc))-datenum(end_date))/366 1];
            rate_shifts_immi = [7/366:7/366:(datenum(sample_cutoff(sc))-datenum(end_date))/366 1];
    
            f = fopen(sprintf('../simulation_results/eir_%d.tsv', rep));
            c=1;
            id = cell(0,0);
            date = cell(0,0);
            date_val = zeros(0,0);
            local_clusters = cell(1000,1);
    
            if strcmp(scheme{v}, 'half')
                dat = textscan(f, '%s\t%s\t%s\t%s\n','HeaderLines',1);fclose(f);
                all_seqs = dat{1};
                subsampled_seqs = randsample(all_seqs, length(all_seqs)/2);
            end 
            f = fopen(sprintf('../simulation_results/eir_%d.tsv', rep));
            while ~feof(f)
                line = strsplit(fgets(f), '\t');
                date_num = datenum(line{2});
                if date_num<=datenum(sample_cutoff(sc))
                    if strcmp(scheme{v}, 'half')
                        test = ismember(line{1}, subsampled_seqs);
                        if test
                            id{c,1} = line{1};
                            date{c,1} = line{2};
                            local_clusters{c,1} = str2double(line{3});
                            
                            %local_clusters{str2double(line{3})} = [local_clusters{str2double(line{3})} ',' id{c,1}];
                            c=c+1;
                        end
                    elseif strcmp(scheme{v}, 'all')
                        id{c,1} = line{1};
                        date{c,1} = line{2};
                        local_clusters{c,1} = str2double(line{3});
                        
                        %local_clusters{str2double(line{3})} = [local_clusters{str2double(line{3})} ',' id{c,1}];
                        c=c+1;
                    end
                end
            end
            local_clusters = local_clusters(~cellfun('isempty',local_clusters));
            local_clusters_numeric = cell2mat(local_clusters);
            local_clusters = unique(local_clusters_numeric);
            % for i = 1:length(local_clusters)
            %     local_clusters{i} = local_clusters{i}(2:end);
            % end
            fclose(f);
    
    
            f = fopen('../../multitree_coalescent/templates/multicoal_template.xml');
            g = fopen(sprintf('../xmls/simmulticoal_%d_%s.xml', rep, scheme{v}), 'w');
            s =  fopen(sprintf('../simulation_results/eir_%d.fasta', rep));
            fgets(s);c=1;seq_id=cell(0,0);
            while ~feof(s)
                line = strsplit(strtrim(fgets(s)));
                fasta(c).Header = line{1};
                seq_names{c} = line{1};
                fasta(c).Sequence = line{2};
                c = c+1;
            end
            fclose(s);
           
            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_data')
                    for lc = 1:length(local_clusters)
                        names = id(ismember(local_clusters_numeric, local_clusters(lc)));
                       % disp(names)
                        fprintf(g, '\t<data id="S%d" spec="Alignment" name="alignment">\n', local_clusters(lc));            
                        for k = 1 :length(names)
                            ind_name = find(ismember(seq_names, names{k}));
                            %disp(ind_name)
                            fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',names{k},names{k}, fasta(ind_name).Sequence);    
                        end
                        fprintf(g, '\t</data>\n');
    
                    end
    
                elseif contains(line, 'insert_tree')
                    for lc = 1:length(local_clusters)
                        fprintf(g, '\t\t<tree id="Tree.t:S%d" spec="beast.evolution.tree.Tree" name="stateNode">\n', local_clusters(lc));                            
                        dates = 'rem';
                        names = id(ismember(local_clusters_numeric, local_clusters(lc)));
                        for k = 1 :length(names)
                            ind = find(ismember(seq_names, names{k}));
                            ind_date = find(ismember(id, names{k}));
                            raw_date = date{ismember(id, names{k})};
                            split_date = strsplit(raw_date,'-');
                            dates = [dates ',' names{k} '=' num2str(str2double(split_date{1}) +(datenum(raw_date, 'yyyy-mm-dd') - datenum(split_date{1}, 'yyyy'))/(datenum(num2str(str2double(split_date{1})+1), 'yyyy') - datenum(split_date{1}, 'yyyy'))) ];
                            %dates = [dates ',' names{k} '=' date{ind_date}];
                        end
                        dates = strrep(dates, 'rem,','');
                        fprintf(g, '\t\t\t<trait id="dateTrait.t:S%d" spec="beast.evolution.tree.TraitSet" dateFormat="decimal" traitname="date" value="%s">\n', local_clusters(lc), dates);
                        fprintf(g, '\t\t\t\t<taxa id="TaxonSet.S%d" spec="TaxonSet">\n', local_clusters(lc));
                        fprintf(g, '\t\t\t\t\t<alignment idref="S%d"/>\n', local_clusters(lc));
                        fprintf(g, '\t\t\t\t</taxa>\n');
                        fprintf(g, '\t\t\t</trait>\n');
                        fprintf(g, '\t\t\t<taxonset idref="TaxonSet.S%d"/>\n', local_clusters(lc));
                        fprintf(g, '\t\t</tree>\n');
                  %      fprintf(g, '\t\t<parameter id="rootHeight:S%s" spec="parameter.RealParameter" lower="0.0" upper="0.25" name="stateNode">0.01</parameter>\n', local_clusters(lc));
                    end
    
                elseif contains(line, 'insert_rate_shifts')
                    for lc = 1:length(local_clusters)
                        fprintf(g, '\t\t<parameter id="rootLength:S%d" name="stateNode" upper="1.0" dimension="1">0.001</parameter>\n', local_clusters(lc));
                    end
                     fprintf(g,'\t\t<parameter id="rateShifts" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts));
                     fprintf(g,'\t\t<parameter id="rateShifts.immi" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts_immi));
    
                elseif contains(line, 'insert_init_tree')
                    for lc = 1:length(local_clusters)
                        fprintf(g, '\t\t<init spec="beast.util.ClusterTree" id="RandomTree.t:S%d" initial="@Tree.t:S%d" clusterType="upgma" taxa="@S%d"/>\n', local_clusters(lc),local_clusters(lc),local_clusters(lc));
    
    
                    end
    
                elseif contains(line, 'insert_priors')
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior2" name="distribution" x="@sigma.immi">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.4" name="distr" meanInRealSpace="true" M="0.5" S="0.25"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="20" S="0.5"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');
              %      fprintf(g,'\t\t\t\t</distribution>\n');
                    fprintf(g,'\t\t\t\t<distribution spec=''nab.skygrid.GrowthRateSmoothingPriorRealParam'' NeLog="@Ne" rateShifts="@rateShifts">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
                    fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-5" sigma="5"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');
                    fprintf(g,'\t\t\t\t<distribution spec=''beast.mascot.skyline.LogSmoothingPrior'' NeLog="@immigrationRate">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.immi"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
                    fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="0" sigma="10"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');
                    fprintf(g,'\t\t\t\t<distribution id="CoalescentConstant.t" spec="nab.multitree.MultiTreeCoalescent" rateIsBackwards="true">\n');
                    fprintf(g,'\t\t\t\t\t<populationModel id="Skygrid" spec="nab.skygrid.Skygrowth" logNe="@Ne" rateShifts="@rateShifts"/>\n');
                    fprintf(g,'\t\t\t\t\t<immigrationRate id="timeVaryingMigrationRates" spec="nab.skygrid.TimeVaryingRates" rate="@immigrationRate" rateShifts="@rateShifts.immi"/>\n');
                    fprintf(g,'\t\t\t\t\t<multiTreeIntervals id="TreeIntervals.t" spec="nab.multitree.MultiTreeIntervals">\n');
                    for lc = 1:length(local_clusters)
                  %      cluster_dates = dat{3}(ismember(dat{2},local_clusters(lc)));
                   %     max_cluster_date = max(datenum(cluster_dates, 'yyyy-mm-dd'));
                      %  mod_cluster_date = datenum(max_cluster_date, 'yyyy-mm-dd');
                    %    offset = ((max_date)-max_cluster_date)/365;
                        fprintf(g, '\t\t\t\t\t\t<tree idref="Tree.t:S%d"/>\n', local_clusters(lc));  
                     %   fprintf(g,'\t\t\t\t\t\t<parameter id="offset:S%s" estimate="true" name="offset">%f</parameter>\n', local_clusters(lc), offset);
                        fprintf(g, '\t\t\t\t\t\t<rootLength idref="rootLength:S%d"/>\n', local_clusters(lc));
                    end
                    
    
                elseif contains(line, 'insert_likelihood')
                    for lc = 1:length(local_clusters)
                        names = id(ismember(local_clusters_numeric, local_clusters(lc)));
                        if length(names)>1
                                fprintf(g, '\t\t\t\t<distribution id="treeLikelihood.S%d" spec="ThreadedTreeLikelihood" data="@S%d" tree="@Tree.t:S%d" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n', local_clusters(lc), local_clusters(lc), local_clusters(lc));                
    
                        end
    
                    end
                    
    
    
    
                elseif contains(line, 'insert_operators')
                     max_size = 0;
                     for a = 1 : length(local_clusters)
                        names = id(ismember(local_clusters_numeric, local_clusters(a)));
                        max_size = max([max_size, length(names)]);
                     end
     
    
                    for lc = 1:length(local_clusters)
                        names = id(ismember(local_clusters_numeric, local_clusters(lc)));
                        if length(names)>1
                            weight = length(names)/max_size;
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:S%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc), 3*weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:S%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc), 3*weight);
                           % fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:S%s" spec="Uniform" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters(lc),local_clusters(lc),30*weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:S%d" spec="SubtreeSlide" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc),45*weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:S%d" spec="Exchange" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc),3*weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:S%d" spec="Exchange" isNarrow="True" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc), 15*weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:S%d" spec="WilsonBalding" tree="@Tree.t:S%d" weight="%f"/>\n',local_clusters(lc),local_clusters(lc), 1*weight);
                        end
                        fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:S%d" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:S%d" weight="0.1"/>\n',local_clusters(lc),local_clusters(lc));
    
                    end
    
                elseif contains(line, 'insert_logs')
                    fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');
                    fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');
    
                     for lc = 1:length(local_clusters)
                        fprintf(g, '\t\t\t<log idref="rootLength:S%d"/>\n', local_clusters(lc));        
                     end
    
                     for lc = 1 : length(local_clusters)
                        fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:S%d"/>\n',lc,local_clusters(lc));
                  %      fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger1:%d" spec="nab.util.MultiTreeStatLogger" heightOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters(lc),local_clusters(lc));
                   %     fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger2:%d" spec="nab.util.MultiTreeStatLogger" originOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters(lc),local_clusters(lc));
                     end
                 elseif contains(line, 'insert_logtree')
    
                    fprintf(g,'\t\t<logger id="typedTreelogger.t:S1" spec="Logger" fileName="$(filebase).trees" logEvery="500000" mode="tree">\n');
                    fprintf(g,'\t\t\t<log id="structuredTreelog.t:S1" spec="nab.multitree.MultiTreeLogger" multiTreeIntervals="@TreeIntervals.t"/>\n');
                    fprintf(g,'\t\t</logger>\n');
                    
    
                else
                    fprintf(g, line);
                end
            end
        end
    end
end
    
   