clc;
clear all;

format long;
format compact;
for problem_size = [10 30 50]
    
    max_nfes = 10000 * problem_size;
    
    rand('seed', sum(100 * clock));
    
    val_2_reach = 10^(-8);
    max_region = 100.0;
    min_region = -100.0;
    lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
    fhd=@cec17_func;
    optimum = 100:100:3000;
    f_name = ['FD-DE_', num2str(problem_size),'_D', '.txt'];
    f_out = fopen(f_name, "wt");
    
    file_name_locMedian = ['FD-DE_',num2str(problem_size), 'D_locMedian.txt'];
    file_locMedian = fopen(file_name_locMedian, "wt");
    for func = 1:30
        outcome = [];
        fprintf(f_out, ['fid:', num2str(func), '\n']);
        disp(['fid:', num2str(func)]);
        
        totalTime = 51;
        for run_id = 1:totalTime
            %%  parameter settings for L-SHADE
            p_best_rate = 0.11;
            arc_rate = 1.4;
            memory_size = 5;
            pop_size = round(25*log(problem_size)*sqrt(problem_size));
            
            max_pop_size = pop_size;
            min_pop_size = 4.0;
            
            %% Initialize the main population
            popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* (repmat(lu(2, :) - lu(1, :), pop_size, 1));
            pop = popold; % the old population becomes the current population
            
            fitness = feval(fhd,pop',func);
            fitness = fitness';
            
            nfes = 0;
            %%%%%%%%%%%%%%%%%%%%%%%% for out
            [bsf_fit_var,gbestid] = min(fitness);
            bsf_solution = pop(gbestid,:);
            memory_old = 0.5 .* ones(memory_size, 1);
            memory_old_cr = 0.5 .* ones(memory_size, 1);
            memory_cr = 0.8 .* ones(memory_size, 1);
            memory_sf = 0.5.*ones(memory_size,1);
            memory_pos = 1;
            
            archive.NP = arc_rate * pop_size; % the maximum size of the archive
            archive.pop = zeros(0, problem_size); % the solutions stored in te archive
            archive.funvalues = zeros(0, 1); % the function value of the archived solutions
            n =   0.6 * pop_size * problem_size;
            counter = zeros(pop_size,1);
           % m = ones(1,2);
            V_lim = 1;
            for i=1:problem_size
                V_lim = V_lim*(abs(lu(1,i)-lu(2,i)));
            end
            V_lim = sqrt(V_lim);
            itter_record_name = ['FD-DE_fid_',num2str(func),'_',num2str(problem_size),'D_',num2str(run_id),'.dat'];
            f_itter_out = fopen(itter_record_name,'wt');
            print_interval = (max_nfes)/20;
            print_num = print_interval;
            fprintf(f_itter_out,'%d\t%.15f\n',1,bsf_fit_var-optimum(func));
             g = 1;
            %% main loop
            while nfes < max_nfes
                pop = popold; % the old population becomes the current population
                [temp_fit, sorted_index] = sort(fitness, 'ascend');
                
                mem_rand_index = ceil(memory_size * rand(pop_size, 1));
                mu_sf = memory_sf(mem_rand_index);
                mu_cr = memory_cr(mem_rand_index);
                
                %% for generating crossover rate
               cr = normrnd(mu_cr, 0.1);
                term_pos = find(mu_cr == -1);
                cr(term_pos) = 0;
                if nfes < 0.4 * max_nfes
                    cr = max(cr,0.6);
                else
                    cr = min(cr, 1);
                    cr = max(cr, 0);
                end
                
                %% for generating scaling factor
                if(nfes < max_nfes * 0.2)
                    sf = sqrt(2).*pi^(-1/3).*(1-mu_sf.^2).*exp(-mu_sf.^2/2) + 0.1 * sin(pi * (rand(pop_size, 1) - 0.8));
                    sf = min(sf, 0.6);
                else
                    sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
                    
                    pos = find(sf <= 0);
                    
                    while ~ isempty(pos)
                        sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                        pos = find(sf <= 0);
                    end
                end
                sf = min(sf, 1);
                
                r0 = [1 : pop_size];
                popAll = [pop; archive.pop];
                [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
                
                pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
                randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
                randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
                pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
                
                vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
                vi = boundConstraint(vi, pop, lu);
                
                mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
                rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; % choose one position where the element of ui doesn't come from the parent
                jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
                ui = vi;
                t1 = rand;
                 rand_mask = rand(pop_size,problem_size);
                bsf_solution_mask = bsf_solution(ones(pop_size, 1), :);
                if t1 > 0.05
                    ui(mask)  = pop(mask);
                else
                    ui(mask) = pop(mask) + rand_mask(mask).*(tpdf(g,1)+ 1).*std(bsf_solution);
                end
                
                
                children_fitness = feval(fhd, ui', func);
                children_fitness = children_fitness';
                
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                for i = 1 : pop_size
                    nfes = nfes + 1;
                    
                    if children_fitness(i) < bsf_fit_var
                        bsf_fit_var = children_fitness(i);
                        bsf_solution = ui(i, :);
                    end
                    if nfes > max_nfes; break; end
                end
                %%%%%%%%%%%%%%%%%%%%%%%% for out
                
                 X = abs(fitness - children_fitness);
                dif = abs((X-mean(X))./(X));
                
                
                %% I == 1: the parent is better; I == 2: the offspring is better
                I = (fitness > children_fitness);
                goodCR = cr(I == 1);
                goodF = sf(I == 1);
                dif_val = dif(I == 1);
                
                %      isempty(popold(I == 1, :))
                archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
                
                [fitness, I] = min([fitness, children_fitness], [], 2);
                
                popold = pop;
                popold(I == 2, :) = ui(I == 2, :);
                
                num_success_params = numel(goodCR);
                
                if num_success_params > 0
                    sum_dif = sum(dif_val);
                    dif_val = dif_val / sum_dif;
                    
                    %% for updating the memory of scaling factor
                    memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
                    memory_sf(memory_pos) = (memory_sf(memory_pos) + memory_old(memory_pos))/2;
                    memory_old(memory_pos) = memory_sf(memory_pos);
                    
                    %% for updating the memory of crossover rate
                    if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
                        memory_cr(memory_pos)  = -1;
                    else
                        memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                    end
                    
                    memory_pos = memory_pos + 1;
                    if memory_pos > memory_size;  memory_pos = 1; end
                end
                C_o = 0;
                for i = 1:pop_size
                    if I(i) == 1
                        counter(i) = counter(i) + 1;
                        C_o = C_o + counter(i);
                    else
                        counter(i) = 0;
                    end
                end
                V_pop = 1;
                for i=1:problem_size
                    V_pop = V_pop.*abs((max(popold(:,i)-min(popold(:,i)))))./2;
                end
                V_pop = sqrt(V_pop);
                nVOL = sqrt(V_pop/V_lim);
                if nVOL < 0.01 
                    for i = 1:pop_size
                        if  C_o >  n 
                            
                            nDim = randi(problem_size);
                            nSeq=randperm(problem_size);
                            j=nSeq(1:nDim);
                            [~,j_size] = size(j);
                            for j1 = 1:j_size
                                j_num = j(1,j1);
                                pop_num = pop(:,j_num);
                                pop_num_1 = pop_num(randperm(numel(pop_num),1));
                                pop(i,j_num)=pop_num_1;
                            end
                            fitness(i) = feval(fhd,pop(i,:)',func);
	            nfes= nfes+1;
                            counter(i)=0;
                        end
                    end
                end
                %% for resizing the population size
                plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
                
                if pop_size > plan_pop_size
                    reduction_ind_num = pop_size - plan_pop_size;
                    if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
                    
                    pop_size = pop_size - reduction_ind_num;
                    for r = 1 : reduction_ind_num
                        [valBest indBest] = sort(fitness, 'ascend');
                        worst_ind = indBest(end);
                        popold(worst_ind,:) = [];
                        pop(worst_ind,:) = [];
                        fitness(worst_ind,:) = [];
	        counter(worst_ind,:) = [];
                    end
                    
                    archive.NP = round(arc_rate * pop_size);
                    
                    if size(archive.pop, 1) > archive.NP
                        rndpos = randperm(size(archive.pop, 1));
                        rndpos = rndpos(1 : archive.NP);
                        archive.pop = archive.pop(rndpos, :);
                    end
                end
                
                if nfes >= print_num
                    fprintf(f_itter_out,'%d\t%.15f\n',nfes,bsf_fit_var-optimum(func));
                    print_num = print_num + print_interval;
                end
                g = g + 1;
            end
            
            bsf_error_val = bsf_fit_var - optimum(func);
            
            
            outcome = [outcome bsf_error_val];
            
            disp(['x=[', num2str(bsf_solution),']=', num2str(bsf_error_val)]);
            fprintf(f_out,'x[%s]=%s\n',num2str(bsf_solution),num2str(bsf_error_val));
            fclose(f_itter_out);
        end %% end 1 run
        fprintf(f_out, ['mean[',num2str(func) , ']=', num2str(mean(outcome)), '\n']);
        fprintf(f_out, ['median[',num2str(func) , ']=', num2str(mean(outcome)), '\n']);
        fprintf(f_out, ['std[',num2str(func) , ']=', num2str(std(outcome)), '\n']);
        disp(['mean[',num2str(func) , ']=', num2str(mean(outcome))]);
        disp(['median[',num2str(func) , ']=', num2str(median(outcome))]);
        disp(['std[',num2str(func) , ']=', num2str(std(outcome))]);
        
        [~,index] = sort(outcome);
        if mod(totalTime,2)==0
            locMedian = [index(totalTime/2),index(totalTime/2+1)];
        else
            locMedian = index(ceil(totalTime/2));
        end
        fprintf(file_locMedian,'fid:%d\n',func);
        fprintf(file_locMedian,'locMedian = %s\n',num2str(locMedian));
    end %% end 1 function run
    fclose(f_out);
    fclose(file_locMedian);
end