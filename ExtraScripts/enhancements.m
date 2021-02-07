figure

temperature_sum = 0;
for i = 1:numSteps
    % Check for scattering
    randNums = rand(numElectrons,1);
    scattered = find(randNums(:,1) < p_scat);
    if(length(scattered)) %if at least one position is beyond region length
        for k = 1:length(scattered)
            vel.x(scattered(k),1) = v_th/sqrt(2) * randn(); % assign new velocity
            vel.y(scattered(k),1) = v_th/sqrt(2) * randn();
        end
    end
    
    % Update x positions
    pos.x(:,2) = pos.x(:,1) + vel.x(:,1) * dt; % (:,1) = old position, get new pos
    
    % Accommodate x=regionLength BC
    pastRightBC = find(pos.x(:,2) > regionLength); % element = 1 if beyond BC
    if(length(pastRightBC)) %if at least one position is beyond region length
        for k = 1:length(pastRightBC)
            pos.x(pastRightBC(k),1) = pos.x(pastRightBC(k),1) - regionLength;
            pos.x(pastRightBC(k),2) = pos.x(pastRightBC(k),2) - regionLength;
        end
    end
    
    % Accommodate x=0 BC
    pastLeftBC = find(pos.x(:,2) < 0); % element = 1 if beyond BC
    if(length(pastLeftBC)) %if at least one position is beyond region length
        for k = 1:length(pastLeftBC)
            pos.x(pastLeftBC(k),1) = pos.x(pastLeftBC(k),1) + regionLength;
            pos.x(pastLeftBC(k),2) = pos.x(pastLeftBC(k),2) + regionLength;
        end
    end
    
    % Update y positions
    pos.y(:,2) = pos.y(:,1) + vel.y(:,1) * dt; % (:,1) = old position, get new pos
    
    % Accommodate y=0 or y=regionWidth BC's
    pastTopOrBottomBC = find((pos.y(:,2) > regionWidth) | (pos.y(:,2) < 0)); % find electrons past top/bottom BC
    if(length(pastTopOrBottomBC)) %if at least one position is beyond region length
        for k = 1:length(pastTopOrBottomBC)
            vel.y(pastTopOrBottomBC(k)) = -1 * vel.y(pastTopOrBottomBC(k)); %angle in = angle out
        end
    end
    
    % Accommodate the boxes' BC's (hard-coded boundaries are easier)
    boxesCheck = find( ((pos.x(:,2) > 80e-9) & (pos.x(:,2) < 120e-9)) & ((pos.y(:,2) < 40e-9) | (pos.y(:,2) > 60e-9)));
    if(length(boxesCheck))
        for k = 1:length(boxesCheck)
            if(pos.y(boxesCheck(k),2) >= 60e-9) % Inside box 1
                if(boxes.box1.type == 'Specular') % Specular Case
                    % Left side
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Bottom
                    else
                        vel.y(boxesCheck(k),1) = -1 * vel.y(boxesCheck(k),1);
                    end
                else % Diffusive case
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        % Only way to hit left wall is with positive x vel
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) > 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        % Must have negative x vel to hit the right side
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) < 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate positive random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Bottom
                    else
                        % Must have positive y vel to hit the bottom side
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.y(boxesCheck(k),1) > 0)
                            % Ensure new y velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    end
                end
            else % Inside box 2
                if(boxes.box1.type == 'Specular') % Specular Case
                    % Left side
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        vel.x(boxesCheck(k),1) = -1 * vel.x(boxesCheck(k),1);
                    % Top
                    else
                        vel.y(boxesCheck(k),1) = -1 * vel.y(boxesCheck(k),1);
                    end
                else % Diffusive case
                    if(pos.x(boxesCheck(k),1) < 80e-9)
                        % Only way to hit left wall is with positive x vel
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) > 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Right side
                    elseif(pos.x(boxesCheck(k),1) > 120e-9)
                        % Must have negative x vel to hit the right side
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.x(boxesCheck(k),1) < 0)
                            % Ensure new x velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate positive random vel
                        end
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    % Bottom
                    else
                        % Must have positive y vel to hit the bottom side
                        vel.y(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                        while(vel.y(boxesCheck(k),1) > 0)
                            % Ensure new y velocity is not into the box
                            vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn(); %Generate negative random vel
                        end
                        vel.x(boxesCheck(k),1) = v_th/sqrt(2) * randn();
                    end
                end
            end
        end
    end
    
    % Re-update positions to ensure new velocities are used for new pos's
    pos.x(:,2) = pos.x(:,1) + vel.x(:,1) * dt; % (:,1) = old position, get new pos
    pos.y(:,2) = pos.y(:,1) + vel.y(:,1) * dt; % (:,1) = old position, get new pos
    
    % Calculate temperature
    % (1/2) m (vel^2)_avg = (1/2)kT * #DoF (2D = 2 DoF
    velocity(:,1) = sqrt(vel.x(:,1).^2 + vel.y(:,1).^2);
    vel_avg = mean(velocity.^2);
    step_temperature = ( (0.5) * (const.m_n) * (vel_avg) ) / (const.k);
    temperature_sum = temperature_sum + step_temperature;
    temperature_avg = temperature_sum/i;
    
    % Plot the electrons
    for k = 1:numToPlot
        subplot(2,1,1);
        plot([pos.x(k,1) pos.x(k,2)], [pos.y(k,1) pos.y(k,2)], 'color', eColours(k,:));
        hold on;
        
        if k == 1
            title('Enhancements')
            xlim([0 regionLength]);
            ylim([0 regionWidth]);
            xlabel('Length [m]');
            ylabel('Width [m]');
            
            rectangle('Position', boxes.box1.coords);
            hold on;
            rectangle('Position', boxes.box2.coords);
            hold on;
        end
    end
    
    %Plot average temperature
    subplot(2,1,2);
    plot(i, step_temperature, 'r.', i, temperature_avg, 'b.');
    title(['Average Temperature: ', num2str(temperature_avg)]);
    hold on;
    if i == 1
        xlim([0 numSteps]);
        xlabel('Step');
        ylabel('Temperature [K]');
    end
    
    % Old positions now equal new positions
    pos.x(:,1) = pos.x(:,2);
    pos.y(:,1) = pos.y(:,2);
    
    pause(0.001)
end

% Electron Distribution Map
figure
hist3([pos.x(:,2), pos.y(:,2)], 'Nbins', [20 10], 'CDataMode','auto')
title('Electron Density')
xlabel('Length [m]'); ylabel('Width [m]');
colorbar
view(2);

% Temperature Distribution Map
temp_dist = zeros(20, 10);
for y = 1:10
    %Bin the y axis
    yhigh = y * 10;
    ylow = yhigh - 10;
    for x = 1:20
        % Bin the x axis
        xhigh = x * 10;
        xlow = xhigh - 10;
        
        % find any electrons within the current bin
        bin = find((pos.x(:,2) > xlow*1e-9) & (pos.x(:,2) < xhigh*1e-9) & (pos.y(:,2) > ylow*1e-9) & (pos.y(:,2) < yhigh*1e-9));
        bin_vel_sum = 0;
        if(~isempty(bin))
            for k = 1:length(bin)
                % Calculate each electron's velocity for averaging
                elec_vel = vel.x(bin(k),1).^2 + vel.y(bin(k),1).^2;
                bin_vel_sum = bin_vel_sum + elec_vel;
            end
            % Calculate temp based on average velocity
            temp_dist(x,y) = ( (0.5) * (const.m_n) * (bin_vel_sum/length(bin)) ) / (const.k);
        else
            % If no electrons in bin, set temp = 0
            temp_dist(x,y) = 0;
        end
        
    end
end
figure
surf(temp_dist');
title('Temperature Distribution');
xlabel('Length [10nm]'); ylabel('Width [10nm]');
colorbar
view(2);


