function pl = plotBloch(cmds)
    clf;
    hold on;
    axis equal;
    axis tight;
    for l=1:length(cmds)
        switch cmds{l}.type
            case 'spline'
                drawSpline(cmds{l});
            case 'equator'
                drawEquator(cmds{l});
            case 'disc'
                drawDisc(cmds{l});
            case 'sphere'
                drawSphere(cmds{l});
            case 'vector'
                drawArrow(cmds{l});
            case 'label'
                drawLabel(cmds{l});
            otherwise
                error('Unknown command ''%s''',cmds{l}.type);
        end
    end
    
    pl = gcf;
end


function drawSphere(cmd)
    cmd=def(cmd,struct('radius',1,'color',[.3 .3 .3 .03]));
    [x y z] = sphere(128);
    r=cmd.radius;
    h = surfl(x*r, y*r, z*r); 
    if length(cmd.color)>3
        color=cmd.color(1:3);
        shade=cmd.color(4);
    else
        color=cmd.color;
        shade=cmd.shade;
    end
    set(h, 'FaceColor',color,'FaceAlpha', shade,'EdgeColor','None');
end

function drawEquator(cmd)
    cmd=def(cmd,struct('radius',1,'dt',0.01,'color',[.3 .3 .3 .03]));
    dt=cmd.dt;
    t=[0:dt:(2*pi) 0];
    r=cmd.radius;
    vecs=r*[cos(t) ; sin(t) ; zeros(1,length(t))];
    c=plot3(vecs(1,:), vecs(2,:), vecs(3,:));
    if length(cmd.color)>3
        color=cmd.color(1:3);
    else
        color=cmd.color;
    end
    set(c,'Color',color);
end

function drawSpline(cmd)
    cmd=def(cmd,struct('radius',1,'dt',0.01,'color',[.3 .3 .3 .03]));
    dt=cmd.dt;
    r=cmd.radius;
    t=[-1*r:dt:r 0];
    vecs=[ zeros(1,length(t)); zeros(1,length(t)) ; t];
    c1=plot3(vecs(1,:), vecs(2,:), vecs(3,:));
    vecs=[ zeros(1,length(t));t; zeros(1,length(t))];
    c2=plot3(vecs(1,:), vecs(2,:), vecs(3,:));
    vecs=[ t; zeros(1,length(t)); zeros(1,length(t))];
    c3=plot3(vecs(1,:), vecs(2,:), vecs(3,:));
    if length(cmd.color)>3
        color=cmd.color(1:3);
    else
        color=cmd.color;
    end
    set(c1,'Color',color);
    set(c2,'Color',color);
    set(c3,'Color',color);
end

function drawArrow(cmd)
    cmd=def(cmd,struct('val',[1 0 0],'radius',1,'color',[.3 .3 .3 .03]));
    r=cmd.radius;
    v = r*1.1*cmd.val;
    if length(cmd.color)>3
        color=cmd.color(1:3);
    else
        color=cmd.color;
    end
    quiver3(0,0,0,v(1),v(2),v(3),'Color',color);
end

function drawLabel(cmd)
    cmd=def(cmd,struct('val',[0,0,1],'label','|0>','radius',1,'color',[.3 .3 .3 .03]));
    r=cmd.radius;
    v=cmd.val;
    text(v(1)*r,v(2)*r,v(3)*r,cmd.label);
end

function s=def(s,def)
    f=fields(def);
    for i=1:length(f)
        l=f{i};
        if ~isfield(s,l) || isempty(s.(l))
            s.(l)=def.(l);
        end
    end
end