clear
close all



Folder = 'TrackData_LKP/Results/';
Folder = 'TrackData_Sochi/Results/';
Folder = 'Straight/Results/';
Folder = 'TrackData_ParkCity/Results/';
Folder = 'TrackData_ParkCity_noentries/Results/';


plot_dims = 0;
plot_spine = 0;
plot_path = 0;
plot_profiles = 0;
plot_surface = 1;

if (plot_dims)
    figure
    dims = load([Folder 'dims.txt']);
    plot(dims(:,1),dims(:,2))
    hold on
    plot(dims(:,1),dims(:,3));
    plot(dims(:,1),0.5*(dims(:,3)-dims(:,2)))

    figure
    plot(dims(:,1),dims(:,4))
    hold on
    plot(dims(:,1),dims(:,5));
end

if (plot_spine +plot_path +plot_surface)
    f2D = figure;
    f3D = figure;
end

if (plot_spine)


    % Plot Track Spine
    track = load([Folder 'spine_spline_pts.txt']);

    yaw = 180*atan2(track(:,6),track(:,5))/pi;
    figure
    plot(track(:,1),yaw);
    xlabel('s')
    ylabel('yaw');

    figure
    plot(track(:,1),track(:,5))
    hold on
    plot(track(:,1),track(:,6))
    plot(track(:,1),track(:,7))
    xlabel('s')
    ylabel('tangent components');

    figure
    plot(track(:,1),track(:,8))
    hold on
    plot(track(:,1),track(:,9))
    plot(track(:,1),track(:,10))
    xlabel('s')
    ylabel('dtangent components');

    figure(f2D)
    plot(track(:,2),track(:,3));
    axis equal
    hold all
    figure(f3D)
    plot3(track(:,2),track(:,3),track(:,4));
    axis equal
    hold all
end

if (plot_path)
    % Plot luge path on track
    path = load([Folder 'withIce_centerline.txt']);
    figure(f2D)
    plot(path(:,2),path(:,3));
    figure(f3D)
    plot3(path(:,2),path(:,3),path(:,4));
end

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

if (plot_surface)
    track_surface = importdata([Folder 'withIce_pts.txt'],' ',1);
    [nUV] = sscanf(track_surface.textdata{1},'%d %d');
    track_grid = reshape(track_surface.data,nUV(2),nUV(1),3);
    figure(f3D)
    surf(track_grid(:,:,1),track_grid(:,:,2),track_grid(:,:,3));
    view([0,0,1])
    axis equal


    npts = size(track_grid,2);
    nsects = 10;
    divs = 1:floor(npts/nsects):npts;
    for i=1:length(divs)-1
        figure
        surf(track_grid(:,divs(i):divs(i+1),1),track_grid(:,divs(i):divs(i+1),2),track_grid(:,divs(i):divs(i+1),3));
        view([0,0,1])
        axis equal
    end

end
