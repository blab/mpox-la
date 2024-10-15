clear;
%start_date = '2020-01-31';
clade = {'all'};
fastafiles = dir('../data/*.fasta');
for i = 1 : length(fastafiles)
    s = fastaread(['../data/' fastafiles(i).name]);
    seq_id=cell(0,0);
    for i = 1 : length(s)
        seq_id{i} = s(1).Header;
    end
    for v = 1: length(clade)
        f = fopen(['../data/' strrep(fastafiles(1).name, 'fasta','tsv')]);
        dat = textscan(f, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n','HeaderLines',1);fclose(f);
        seq_names = cell(0,0);
        for j = 1 : length(s)
            seq_names{j,1} = s(j).Header;
        end


        all_clusters = unique(dat{2});              
    
    
        
        
        
         
         % for a = length(dat{1}):-1:1
         %     clades = dat{5};
         %     u_clade = unique(clades);
         %     if strcmp(clade{v}, 'Alpha') || strcmp(clade{v}, 'Omicron') || strcmp(clade{v}, 'Delta') 
         %         if sum(ismember(dat{5}{a},clade(v))) == 0
         %             dat{1}(a) = [];
         %             dat{2}(a) = [];
         %             dat{3}(a) = [];
         %             dat{4}(a) = [];
         %             dat{5}(a) = [];
         %             dat{6}(a) = [];
         %         end
         %     end
         % end
         % 
                    
        
        % if strcmp(clade{v}, 'subsampled')
        %     subsample_n = 1500;
        % 
        %     date_list = year(datetime(dat{1,4}))*100+ week(datetime(dat{1,4}));
        % 
        %     u_date = unique(date_list);
        %     tmp_names = dat{1};
        %     sub_names = cell(0,0);
        %     counter = 1;
        %     while counter <= subsample_n
        %         samp_date = randsample(u_date,1);
        %         names = tmp_names(ismember(date_list,samp_date));
        %         if ~isempty(names) 
        %             rand_name = randsample(names, 1);
        %             if ~isempty(rand_name)
        %                 sub_names{counter,1} = rand_name;
        %                 ind_name = find(ismember(tmp_names, sub_names{counter,1}));
        %                 tmp_names(ind_name) = [];
        %                 date_list(ind_name) = [];
        %                 counter = counter + 1 ;
        %             end
        %         else 
        % 
        %         end
        %     end 

            % flat_sub_names = [sub_names{:}];
            % for a = length(dat{1}):-1:1
            %     if sum(ismember(dat{1}{a}, flat_sub_names)) == 0
            %         dat{1}(a) = [];
            %         dat{2}(a) = [];
            %         dat{3}(a) = [];
            %         dat{4}(a) = [];
            %         dat{5}(a) = [];
            %         dat{6}(a) = [];
            %     end 
            % end
        end
        

        all_dates = dat{6};
        max_date = max(datenum(all_dates, 'yyyy-mm-dd'));
        min_date = min(datenum(all_dates, 'yyyy-mm-dd'));
        adjusted_start_date = datetime(min_date,'ConvertFrom','datenum') - calmonths(2);
        rate_shift = [7/366:7/366:(max_date-datenum(adjusted_start_date))/366];
        rate_shift_immi = [7/366:7/366:(max_date-datenum(adjusted_start_date))/366];

        local_clusters = unique(dat{2});

        f = fopen('../templates/multicoal_template_cases.xml');
        g = fopen(['../xmls/multicoal_updated_case_prior_' strrep(fastafiles(1).name, 'fasta','xml')],'w');
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_data')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t<data id="S%s" spec="Alignment" name="alignment">\n', local_clusters{lc});                
                    names = dat{1}(ismember(dat{2},local_clusters{lc}));
                        for k = 1 :length(names)
                            ind_name = find(ismember(seq_names, names{k}));
                            fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',...
                                s(ind_name).Header, s(ind_name).Header, s(ind_name).Sequence);
                        end
                        fprintf(g, '\t</data>\n');

                end

            elseif contains(line, 'insert_tree')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t<tree id="Tree.t:S%s" spec="beast.evolution.tree.Tree" name="stateNode">\n', local_clusters{lc});                            
                    dates = 'rem';
                    names = dat{1}(ismember(dat{2},local_clusters{lc}));
                    for k = 1 :length(names)
                        raw_date = dat{6}{ismember(dat{1},names{k})};
                        split_date = strsplit(raw_date,'-');
                        dates = [dates ',' names{k} '=' num2str(str2double(split_date{1}) +(datenum(raw_date, 'yyyy-mm-dd') - datenum(split_date{1}, 'yyyy'))/(datenum(num2str(str2double(split_date{1})+1), 'yyyy') - datenum(split_date{1}, 'yyyy'))) ];
                    end
                    dates = strrep(dates, 'rem,','');
                    fprintf(g, '\t\t\t<trait id="dateTrait.t:S%s" spec="beast.evolution.tree.TraitSet" dateFormat="decimal" traitname="date" value="%s">\n', local_clusters{lc}, dates);
                    fprintf(g, '\t\t\t\t<taxa id="TaxonSet.S%s" spec="TaxonSet">\n', local_clusters{lc});
                    fprintf(g, '\t\t\t\t\t<alignment idref="S%s"/>\n', local_clusters{lc});
                    fprintf(g, '\t\t\t\t</taxa>\n');
                    fprintf(g, '\t\t\t</trait>\n');
                    fprintf(g, '\t\t\t<taxonset idref="TaxonSet.S%s"/>\n', local_clusters{lc});
                    fprintf(g, '\t\t</tree>\n');
              %      fprintf(g, '\t\t<parameter id="rootHeight:S%s" spec="parameter.RealParameter" lower="0.0" upper="0.25" name="stateNode">0.01</parameter>\n', local_clusters{lc});
                end

            elseif contains(line, 'insert_rate_shifts')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t<parameter id="rootLength:S%s" name="stateNode" upper="1.0" dimension="1">0.001</parameter>\n', local_clusters{lc});
                end
                 fprintf(g,'\t\t<parameter id="rateShifts" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shift));
                 fprintf(g,'\t\t<parameter id="rateShifts.immi" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shift_immi));

            elseif contains(line, 'insert_init_tree')
                for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t<init spec="beast.util.ClusterTree" id="RandomTree.t:S%s" initial="@Tree.t:S%s" clusterType="upgma" taxa="@S%s"/>\n', local_clusters{lc},local_clusters{lc},local_clusters{lc});


                end

            elseif contains(line, 'insert_priors')
                fprintf(g,'\t\t\t\t<prior id="Sigmaprior2" name="distribution" x="@sigma.immi">\n');
                fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.4" name="distr" meanInRealSpace="true" M="0.5" S="0.25"/>\n');
                fprintf(g,'\t\t\t\t</prior>\n');
                fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="20" S="0.5"/>\n');
                fprintf(g,'\t\t\t\t</prior>\n');
          %      fprintf(g,'\t\t\t\t</distribution>\n');
                fprintf(g,'\t\t\t\t<prior id="NeCasePriorErrors" name = "distribution">\n');
                fprintf(g,'\t\t\t\t\t<x id="ErrorTerms" spec="nab.multitree.NeCasePriorErrorDifference" arg="@Ne" logCases="@cases" overallNeScaler="@scaler"/>\n');
                fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                fprintf(g,'\t\t\t\t\t</distr>\n');
                fprintf(g,'\t\t\t\t</prior>\n');
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
              %      cluster_dates = dat{3}(ismember(dat{2},local_clusters{lc}));
               %     max_cluster_date = max(datenum(cluster_dates, 'yyyy-mm-dd'));
                  %  mod_cluster_date = datenum(max_cluster_date, 'yyyy-mm-dd');
                %    offset = ((max_date)-max_cluster_date)/365;
                    fprintf(g, '\t\t\t\t\t\t<tree idref="Tree.t:S%s"/>\n', local_clusters{lc});  
                 %   fprintf(g,'\t\t\t\t\t\t<parameter id="offset:S%s" estimate="true" name="offset">%f</parameter>\n', local_clusters{lc}, offset);
                    fprintf(g, '\t\t\t\t\t\t<rootLength idref="rootLength:S%s"/>\n', local_clusters{lc});
                end
                

            elseif contains(line, 'insert_likelihood')
                for lc = 1:length(local_clusters)
                    names = dat{1}(ismember(dat{2},local_clusters{lc}));
                    if length(names)>1
                            fprintf(g, '\t\t\t\t<distribution id="treeLikelihood.S%s" spec="ThreadedTreeLikelihood" data="@S%s" tree="@Tree.t:S%s" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n', local_clusters{lc}, local_clusters{lc}, local_clusters{lc});                

                    end

                end
                



            elseif contains(line, 'insert_operators')

                for lc = 1:length(local_clusters)
                    names = dat{1}(ismember(dat{2},local_clusters{lc}));
                    if length(names)>1
                        weight = length(names)/length(dat{2})*5; 
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:S%s" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc}, 3*weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:S%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc}, 3*weight);
                       % fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:S%s" spec="Uniform" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc},30*weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:S%s" spec="SubtreeSlide" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc},45*weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:S%s" spec="Exchange" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc},3*weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:S%s" spec="Exchange" isNarrow="True" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc}, 15*weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:S%s" spec="WilsonBalding" tree="@Tree.t:S%s" weight="%f"/>\n',local_clusters{lc},local_clusters{lc}, 1*weight);
                    end
                    fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:S%s" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:S%s" weight="0.1"/>\n',local_clusters{lc},local_clusters{lc});

                end

            elseif contains(line, 'insert_logs')
                fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');
                fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');
                fprintf(g,'\t\t\t<log idref="ErrorTerms"/>\n');

                 for lc = 1:length(local_clusters)
                    fprintf(g, '\t\t\t<log idref="rootLength:S%s"/>\n', local_clusters{lc});        
                 end

                 for lc = 1 : length(local_clusters)
                    fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:S%s"/>\n',lc,local_clusters{lc});
              %      fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger1:%d" spec="nab.util.MultiTreeStatLogger" heightOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters{lc},local_clusters{lc});
               %     fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger2:%d" spec="nab.util.MultiTreeStatLogger" originOnly="true" tree="@Tree.t:S%s"  rootLength="@rootLength:S%s" />\n',lc,local_clusters{lc},local_clusters{lc});
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
