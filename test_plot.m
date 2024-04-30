clear
close all

f2D = figure;
f3D = figure;

Folder = 'TrackData_LKP/Results/';
Folder = 'TrackData_Sochi/Results/';

% Plot Track Spine
track = load([Folder 'track_test_pts.txt']);
figure(f2D)
plot(track(:,1),track(:,3));
axis equal
hold all
figure(f3D)
plot3(track(:,1),track(:,2),track(:,3));
axis equal
hold all

% Plot luge path on track
path = load([Folder 'path.txt']);
figure(f2D)
plot(path(:,2),path(:,4));
figure(f3D)
plot3(path(:,2),path(:,3),path(:,4));

plot_profiles = 0;
if (plot_profiles) 
    % Load profiles
    nfiles = 750;
    for cnt=1:nfiles
        profiles{cnt} = load([Folder +'profile' num2str(cnt) '.pts']);
    end

    figure(f2D)
    for cnt=1:nfiles
        plot(profiles{cnt}(:,1),profiles{cnt}(:,3));
    end

    figure(f3D)
    for cnt=1:nfiles
        plot3(profiles{cnt}(:,1),profiles{cnt}(:,2),profiles{cnt}(:,3));
    end
end

plot_surface = 1;
if (plot_surface)
    track_surface = importdata([Folder 'track_surface_pts.txt'],' ',1);
    [nUV] = sscanf(track_surface.textdata{1},'%d %d');
    track_grid = reshape(track_surface.data,nUV(2),nUV(1),3)
    figure(f3D)
    surf(track_grid(:,:,1),track_grid(:,:,2),track_grid(:,:,3));
end
