%% 2D CPU/GPU Die Cooling Simulation: Explicit & Implicit Transients + Steady State
% Models a chip package cross-section undergoing transient thermal analysis.
%
% Geometry (z = 0 is bottom of silicon die, z = Lz is top of heatspreader):
%   Layer 1 (bottom): Silicon die         with localized heat sources (CPU cores)
%   Layer 2 (middle): TIM                 thermal interface material
%   Layer 3 (top):    Copper heatspreader convection to forced air on top surface
%
% Boundary Conditions:
%   Top surface    (z = Lz): Convection to forced air, h = 2000 W/m^2-K, T_inf = 300 K
%   Left surface   (x = 0):  Convection to air, same h and T_inf
%   Right surface  (x = Lx): Convection to air, same h and T_inf
%   Bottom surface (z = 0):  Adiabatic (insulated PCB side)
%
% Solvers:
%   1. Transient Explicit  (forward-time, central-space, Table 5.3)
%   2. Transient Implicit  (backward-time, central-space, Gauss-Seidel per step)
%   3. Steady State        (Table 4.2 stencils, Gauss-Seidel)
%
% Outputs:
%   - 8 transient snapshots (explicit)
%   - MP4 animation: cpu_cooling_explicit.mp4
%   - Final 3-panel comparison: explicit | steady-state | implicit
%   - Cooldown phase: cores off, chip returns to ambient
%   - MP4 animation: cpu_cooldown_explicit.mp4
%   - Peak temperature vs time plot

clear; clc; close all;

%% =========================================================================
%  GEOMETRY AND GRID
% =========================================================================
Lx = 0.020;   % m, die width  (20 mm)
Lz = 0.006;   % m, total stack height (6 mm)

t_die = 0.003;   % m, silicon die (3 mm)
t_tim = 0.0005;  % m, TIM        (0.5 mm)
t_hs  = 0.0025;  % m, Cu spreader (2.5 mm)

nodes_per_meter = 500;
Nx = round(Lx * nodes_per_meter) + 1;
Nz = round(Lz * nodes_per_meter) + 1;
dx = Lx / (Nx - 1);
dz = Lz / (Nz - 1);

if abs(dx - dz) > 1e-7
    error('Non-square grid: dx=%.4g, dz=%.4g.', dx, dz);
end

x = linspace(0, Lx, Nx);
z = linspace(0, Lz, Nz);
fprintf('Grid: %d x %d nodes | dx = dz = %.4g m\n', Nx, Nz, dx); %[output:5e659e70]

%% =========================================================================
%  MATERIAL PROPERTIES
% =========================================================================
k_si = 130;  rho_si = 2330; cp_si = 700;
k_tim_v = 5.0; rho_tim = 2500; cp_tim = 800;
k_cu = 400;  rho_cu = 8960; cp_cu = 385;

K   = zeros(Nz, Nx);
RHO = zeros(Nz, Nx);
CP  = zeros(Nz, Nx);

for i = 1:Nz
    if z(i) <= t_die
        K(i,:) = k_si;  RHO(i,:) = rho_si;  CP(i,:) = cp_si;
    elseif z(i) <= t_die + t_tim
        K(i,:) = k_tim_v; RHO(i,:) = rho_tim; CP(i,:) = cp_tim;
    else
        K(i,:) = k_cu;  RHO(i,:) = rho_cu;  CP(i,:) = cp_cu;
    end
end

ALPHA = K ./ (RHO .* CP);

%% =========================================================================
%  HEAT SOURCE MAP  (3 CPU cores in bottom of silicon die)
% =========================================================================
Q = zeros(Nz, Nx);
q_core = 5e8;
core_h = t_die * 0.4;
core_w = Lx * 0.12;
core_x_centers = [Lx*0.20, Lx*0.50, Lx*0.80];

[Xg, Zg] = meshgrid(x, z);
for cidx = 1:3
    cx = core_x_centers(cidx);
    Q( Zg <= core_h & abs(Xg - cx) <= core_w/2 ) = q_core;
end

%% =========================================================================
%  BOUNDARY CONDITIONS & TIME STEP
% =========================================================================
T_init = 300;
T_inf  = 300;
h_conv = 2000;

k_min_global = min(K(:));
Bi_global    = h_conv * dx / k_min_global;
alpha_max    = max(ALPHA(:));

Fo_lim = min([1/4, 1/(2*(2+Bi_global)), 1/(4*(1+Bi_global))]);
Fo     = 0.45 * Fo_lim;
dt     = Fo * dx^2 / alpha_max;

fprintf('Fo = %.4g (lim %.4g) | dt = %.4g s | Bi = %.4g\n', Fo, Fo_lim, dt, Bi_global); %[output:8e3bd423]

FO = ALPHA * dt / dx^2;
BI = h_conv * dx ./ K;
QS = Q .* dx^2 ./ K;

%% =========================================================================
%  TRANSIENT EXPLICIT SOLVER  (fully vectorized)
% =========================================================================
T = T_init * ones(Nz, Nx);

tol_abs    = 1e-3;
hold_steps = 10;
max_steps  = 3e6;
pass_ct    = 0;

t_char     = (min(Lx,Lz))^2 / alpha_max;
snap_times = t_char * [1e-4 5e-4 2e-3 8e-3 2e-2 6e-2 1.5e-1];
snap_k = 1; num_snaps = 8; time = 0;

figSnaps = figure('Name','CPU Die Cooling - Explicit Snapshots','Color','w'); %[output:567f3300]
colormap('hot'); %[output:567f3300]

videoName   = 'cpu_cooling_explicit.mp4';
frame_dt    = max(dt, t_char / 400);
nextFrame_t = 0;
vw = VideoWriter(videoName, 'MPEG-4');
vw.FrameRate = 30;
open(vw);

figAnim = figure('Visible','off','Color','w','Position',[100 100 640 400]);
axA  = axes('Parent', figAnim);
hImg = imagesc(axA, x*1e3, z*1e3, T);
set(axA, 'YDir','normal');
colormap(axA,'hot');
cb = colorbar(axA); ylabel(cb,'Temperature (K)');
xlabel(axA,'x (mm)'); ylabel(axA,'z (mm)');
title(axA, sprintf('t = %.4g s', time)); axis(axA,'tight');
hold(axA,'on');
yline(axA, t_die*1e3,         'c--','LineWidth',1.2);
yline(axA, (t_die+t_tim)*1e3, 'b--','LineWidth',1.2);
hold(axA,'off');
drawnow('expose');
writeVideo(vw, getframe(figAnim));
nextFrame_t = frame_dt;

fprintf('\nRunning explicit solver (vectorized)...\n'); %[output:51d49ccd]

for step = 1:max_steps %[output:group:9b10ea1a]

    Tp = T;

    % Interior nodes (Case 1)
    fo = FO(2:Nz-1, 2:Nx-1);
    qs = QS(2:Nz-1, 2:Nx-1);
    Tp(2:Nz-1, 2:Nx-1) = fo.*(T(3:Nz,2:Nx-1) + T(1:Nz-2,2:Nx-1) + ...
                               T(2:Nz-1,3:Nx) + T(2:Nz-1,1:Nx-2)) ...
                        + (1 - 4*fo).*T(2:Nz-1,2:Nx-1) + fo.*qs;

    % Top edge (Case 3)
    i = Nz;
    fo_t = FO(i,2:Nx-1); bi_t = BI(i,2:Nx-1); qs_t = QS(i,2:Nx-1);
    Tp(i,2:Nx-1) = fo_t.*(2*T(i-1,2:Nx-1) + T(i,3:Nx) + T(i,1:Nx-2) + 2*bi_t*T_inf) ...
                 + (1 - 4*fo_t - 2*bi_t.*fo_t).*T(i,2:Nx-1) + fo_t.*qs_t;

    % Bottom edge (adiabatic ghost node)
    i = 1;
    fo_b = FO(i,2:Nx-1); qs_b = QS(i,2:Nx-1);
    Tp(i,2:Nx-1) = fo_b.*(2*T(i+1,2:Nx-1) + T(i,3:Nx) + T(i,1:Nx-2)) ...
                 + (1 - 4*fo_b).*T(i,2:Nx-1) + fo_b.*qs_b;

    % Left edge (Case 3)
    j = 1;
    fo_l = FO(2:Nz-1,j); bi_l = BI(2:Nz-1,j); qs_l = QS(2:Nz-1,j);
    Tp(2:Nz-1,j) = fo_l.*(2*T(2:Nz-1,j+1) + T(3:Nz,j) + T(1:Nz-2,j) + 2*bi_l*T_inf) ...
                 + (1 - 4*fo_l - 2*bi_l.*fo_l).*T(2:Nz-1,j) + fo_l.*qs_l;

    % Right edge (Case 3)
    j = Nx;
    fo_r = FO(2:Nz-1,j); bi_r = BI(2:Nz-1,j); qs_r = QS(2:Nz-1,j);
    Tp(2:Nz-1,j) = fo_r.*(2*T(2:Nz-1,j-1) + T(3:Nz,j) + T(1:Nz-2,j) + 2*bi_r*T_inf) ...
                 + (1 - 4*fo_r - 2*bi_r.*fo_r).*T(2:Nz-1,j) + fo_r.*qs_r;

    % Corners (Case 4)
    i=Nz; j=1;   fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
    Tp(i,j) = 2*fo_c*(T(i-1,j)+T(i,j+1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T(i,j) + 2*fo_c*qs_c;
    i=Nz; j=Nx;  fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
    Tp(i,j) = 2*fo_c*(T(i-1,j)+T(i,j-1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T(i,j) + 2*fo_c*qs_c;
    i=1;  j=1;   fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
    Tp(i,j) = 2*fo_c*(T(i+1,j)+T(i,j+1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T(i,j) + 2*fo_c*qs_c;
    i=1;  j=Nx;  fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
    Tp(i,j) = 2*fo_c*(T(i+1,j)+T(i,j-1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T(i,j) + 2*fo_c*qs_c;

    dT_max = max(abs(Tp(:) - T(:)));
    T    = Tp;
    time = time + dt;

    % Animation frame
    if time >= nextFrame_t
        set(hImg, 'CData', T);
        title(axA, sprintf('t = %.4g s', time));
        drawnow('expose');
        writeVideo(vw, getframe(figAnim));
        nextFrame_t = nextFrame_t + frame_dt;
    end

    % Snapshots 1..7
    figure(figSnaps);
    while (snap_k <= num_snaps-1) && (time >= snap_times(snap_k))
        subplot(2,4,snap_k); %[output:567f3300]
        imagesc(x*1e3, z*1e3, T); set(gca,'YDir','normal');
        colormap('hot'); colorbar;
        title(sprintf('t = %.3g s', time));
        xlabel('x (mm)'); ylabel('z (mm)');
        hold on;
        yline(t_die*1e3,         'c--','LineWidth',1);
        yline((t_die+t_tim)*1e3, 'b--','LineWidth',1);
        hold off;
        snap_k = snap_k + 1;
    end

    pass_ct = pass_ct + (dT_max < tol_abs);
    if dT_max >= tol_abs, pass_ct = 0; end
    if pass_ct >= hold_steps
        fprintf('Explicit converged: t = %.4f s, %d steps, dT_max = %.3g K\n', time, step, dT_max);
        break
    end
end %[output:group:9b10ea1a]

% 8th snapshot
T_jmax = max(T(:));
T_jmin = min(T(:));

figure(figSnaps);
subplot(2,4,num_snaps);
imagesc(x*1e3, z*1e3, T); set(gca,'YDir','normal');
colormap('hot'); colorbar;
title(sprintf('Final (steady)\nt = %.2f s', time));
xlabel('x (mm)'); ylabel('z (mm)');
hold on;
yline(t_die*1e3,         'c--','LineWidth',1);
yline((t_die+t_tim)*1e3, 'b--','LineWidth',1);
hold off;

set(hImg, 'CData', T);
title(axA, sprintf('Final (steady) t = %.3f s', time));
drawnow('expose'); writeVideo(vw, getframe(figAnim));
close(vw); close(figAnim);
fprintf('Saved animation: %s\n', videoName);

T_explicit_final = T;

%% =========================================================================
%  TRANSIENT IMPLICIT SOLVER (vectorized Gauss-Seidel per time step)
% =========================================================================
fprintf('\nRunning implicit solver (vectorized)...\n');

Tn  = T_init * ones(Nz, Nx);
Tnp = Tn;
tol_inner = 1e-6;
max_inner  = 200;
time_imp   = 0;
step_imp   = 0;
pass_ct    = 0;

while step_imp < max_steps
    step_imp = step_imp + 1;
    Tnp = Tn;

    for it = 1:max_inner
        Told_all = Tnp;

        fo = FO(2:Nz-1,2:Nx-1); qs = QS(2:Nz-1,2:Nx-1);
        Tnp(2:Nz-1,2:Nx-1) = (Tn(2:Nz-1,2:Nx-1) ...
            + fo.*(Tnp(3:Nz,2:Nx-1)+Tnp(1:Nz-2,2:Nx-1)+Tnp(2:Nz-1,3:Nx)+Tnp(2:Nz-1,1:Nx-2)) ...
            + fo.*qs) ./ (1 + 4*fo);

        i = Nz;
        fo_t=FO(i,2:Nx-1); bi_t=BI(i,2:Nx-1); qs_t=QS(i,2:Nx-1);
        Tnp(i,2:Nx-1) = (Tn(i,2:Nx-1) ...
            + fo_t.*(2*Tnp(i-1,2:Nx-1)+Tnp(i,3:Nx)+Tnp(i,1:Nx-2)) ...
            + 2*bi_t.*fo_t*T_inf + fo_t.*qs_t) ./ (1 + 2*fo_t.*(2+bi_t));

        i = 1;
        fo_b=FO(i,2:Nx-1); qs_b=QS(i,2:Nx-1);
        Tnp(i,2:Nx-1) = (Tn(i,2:Nx-1) ...
            + fo_b.*(2*Tnp(i+1,2:Nx-1)+Tnp(i,3:Nx)+Tnp(i,1:Nx-2)) ...
            + fo_b.*qs_b) ./ (1 + 4*fo_b);

        j = 1;
        fo_l=FO(2:Nz-1,j); bi_l=BI(2:Nz-1,j); qs_l=QS(2:Nz-1,j);
        Tnp(2:Nz-1,j) = (Tn(2:Nz-1,j) ...
            + fo_l.*(2*Tnp(2:Nz-1,j+1)+Tnp(3:Nz,j)+Tnp(1:Nz-2,j)) ...
            + 2*bi_l.*fo_l*T_inf + fo_l.*qs_l) ./ (1 + 2*fo_l.*(2+bi_l));

        j = Nx;
        fo_r=FO(2:Nz-1,j); bi_r=BI(2:Nz-1,j); qs_r=QS(2:Nz-1,j);
        Tnp(2:Nz-1,j) = (Tn(2:Nz-1,j) ...
            + fo_r.*(2*Tnp(2:Nz-1,j-1)+Tnp(3:Nz,j)+Tnp(1:Nz-2,j)) ...
            + 2*bi_r.*fo_r*T_inf + fo_r.*qs_r) ./ (1 + 2*fo_r.*(2+bi_r));

        i=Nz; j=1;  fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
        Tnp(i,j)=(Tn(i,j)+2*fo_c*(Tnp(i-1,j)+Tnp(i,j+1))+4*bi_c*fo_c*T_inf+2*fo_c*qs_c)/(1+4*fo_c*(1+bi_c));
        i=Nz; j=Nx; fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
        Tnp(i,j)=(Tn(i,j)+2*fo_c*(Tnp(i-1,j)+Tnp(i,j-1))+4*bi_c*fo_c*T_inf+2*fo_c*qs_c)/(1+4*fo_c*(1+bi_c));
        i=1;  j=1;  fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
        Tnp(i,j)=(Tn(i,j)+2*fo_c*(Tnp(i+1,j)+Tnp(i,j+1))+4*bi_c*fo_c*T_inf+2*fo_c*qs_c)/(1+4*fo_c*(1+bi_c));
        i=1;  j=Nx; fo_c=FO(i,j); bi_c=BI(i,j); qs_c=QS(i,j);
        Tnp(i,j)=(Tn(i,j)+2*fo_c*(Tnp(i+1,j)+Tnp(i,j-1))+4*bi_c*fo_c*T_inf+2*fo_c*qs_c)/(1+4*fo_c*(1+bi_c));

        if max(abs(Tnp(:)-Told_all(:))) < tol_inner, break; end
    end

    dT_max   = max(abs(Tnp(:)-Tn(:)));
    Tn       = Tnp;
    time_imp = time_imp + dt;

    pass_ct = pass_ct + (dT_max < tol_abs);
    if dT_max >= tol_abs, pass_ct = 0; end
    if pass_ct >= hold_steps
        fprintf('Implicit converged: t = %.4f s, %d steps, dT_max = %.3g K\n', time_imp, step_imp, dT_max);
        break
    end
end

T_implicit_final = Tn;

%% =========================================================================
%  STEADY-STATE SOLVER (Table 4.2, Gauss-Seidel)
% =========================================================================
fprintf('\nRunning steady-state solver...\n');

Tss = T_init * ones(Nz, Nx);
tol_ss = 1e-3; max_it_ss = 3e6;
res = Inf; it = 0;

while res > tol_ss && it < max_it_ss
    res = 0; it = it + 1;
    for i = 1:Nz
        for j = 1:Nx
            Told   = Tss(i,j);
            bi_loc = BI(i,j);
            qs_loc = QS(i,j);
            if i == 1
                if     j == 1,  Tnew=(Tss(i+1,j)+Tss(i,j+1)+2*bi_loc*T_inf+qs_loc)/(2+2*bi_loc);
                elseif j == Nx, Tnew=(Tss(i+1,j)+Tss(i,j-1)+2*bi_loc*T_inf+qs_loc)/(2+2*bi_loc);
                else,           Tnew=(2*Tss(i+1,j)+Tss(i,j+1)+Tss(i,j-1)+qs_loc)/4;
                end
            elseif i == Nz
                if     j == 1,  Tnew=(Tss(i-1,j)+Tss(i,j+1)+2*bi_loc*T_inf+qs_loc)/(2+2*bi_loc);
                elseif j == Nx, Tnew=(Tss(i-1,j)+Tss(i,j-1)+2*bi_loc*T_inf+qs_loc)/(2+2*bi_loc);
                else,           Tnew=(2*Tss(i-1,j)+Tss(i,j+1)+Tss(i,j-1)+2*bi_loc*T_inf+qs_loc)/(4+2*bi_loc);
                end
            elseif j == 1
                Tnew=(Tss(i-1,j)+Tss(i+1,j)+2*Tss(i,j+1)+2*bi_loc*T_inf+qs_loc)/(4+2*bi_loc);
            elseif j == Nx
                Tnew=(Tss(i-1,j)+Tss(i+1,j)+2*Tss(i,j-1)+2*bi_loc*T_inf+qs_loc)/(4+2*bi_loc);
            else
                Tnew=(Tss(i+1,j)+Tss(i-1,j)+Tss(i,j+1)+Tss(i,j-1)+qs_loc)/4;
            end
            Tss(i,j) = Tnew;
            res = max(res, abs(Tnew-Told));
        end
    end
end
fprintf('Steady-state converged: %d iterations, max change = %.3g K\n', it, res);

%% =========================================================================
%  FINAL COMPARISON FIGURE
% =========================================================================
T_all  = {T_explicit_final, Tss, T_implicit_final};
titles = {'Explicit Transient: Final','Steady State: Final','Implicit Transient: Final'};

figure('Name','Final Comparison','Color','w'); colormap('hot');
for p = 1:3
    subplot(1,3,p);
    imagesc(x*1e3, z*1e3, T_all{p}); set(gca,'YDir','normal');
    set(gca, 'CLim', [T_jmin, T_jmax]);
    colorbar; title(titles{p}); xlabel('x (mm)'); ylabel('z (mm)');
    hold on;
    yline(t_die*1e3,         'c--','LineWidth',1.2,'LabelHorizontalAlignment','left');
    yline((t_die+t_tim)*1e3, 'b--','LineWidth',1.2,'LabelHorizontalAlignment','left');
    hold off;
end

%% =========================================================================
%  KEY RESULTS
% =========================================================================
fprintf('\n--- Key Results ---\n');
fprintf('Peak junction temp : %.2f K  (%.2f deg C)\n', T_jmax, T_jmax-273.15);
fprintf('Min temperature    : %.2f K  (%.2f deg C)\n', T_jmin, T_jmin-273.15);
fprintf('Max temp rise      : %.2f K\n', T_jmax - T_init);
Q_total = sum(Q(:)) * dx^2 * dx;
if Q_total > 0
    fprintf('Thermal resistance : %.4f K/W\n', (T_jmax-T_inf)/Q_total);
end

%% =========================================================================
%  COOLDOWN PHASE (cores off — explicit, vectorized)
% =========================================================================
fprintf('\nRunning cooldown phase (cores off)...\n');

T_cd    = T_explicit_final;
pass_ct = 0;
time_cd = 0;
step_cd = 0;

t_cd_char     = t_char;
snap_times_cd = t_cd_char * [1e-3 5e-3 2e-2 8e-2 2e-1 6e-1 1.5];
snap_k_cd     = 1;

figCool = figure('Name','Cooldown Snapshots (Cores Off)','Color','w');

videoName_cd = 'cpu_cooldown_explicit.mp4';
frame_dt_cd  = max(dt, t_cd_char / 400);
nextFrame_cd = 0;
vw_cd = VideoWriter(videoName_cd, 'MPEG-4');
vw_cd.FrameRate = 30;
open(vw_cd);

figAnimCd = figure('Visible','off','Color','w','Position',[100 100 640 400]);
axC    = axes('Parent', figAnimCd);
hImgCd = imagesc(axC, x*1e3, z*1e3, T_cd);
set(axC, 'YDir','normal', 'CLim', [T_inf, T_jmax]);
colormap(axC,'hot');
cb2 = colorbar(axC); ylabel(cb2,'Temperature (K)');
xlabel(axC,'x (mm)'); ylabel(axC,'z (mm)');
title(axC,'Cooldown: t = 0 s (cores off)'); axis(axC,'tight');
hold(axC,'on');
yline(axC, t_die*1e3,         'c--','LineWidth',1.2);
yline(axC, (t_die+t_tim)*1e3, 'b--','LineWidth',1.2);
hold(axC,'off');
drawnow('expose');
writeVideo(vw_cd, getframe(figAnimCd));
nextFrame_cd = frame_dt_cd;

T_peak_history = max(T_cd(:));
T_time_history = 0;

for step_cd = 1:max_steps

    Tp = T_cd;

    fo = FO(2:Nz-1, 2:Nx-1);
    Tp(2:Nz-1, 2:Nx-1) = fo.*(T_cd(3:Nz,2:Nx-1) + T_cd(1:Nz-2,2:Nx-1) + ...
                               T_cd(2:Nz-1,3:Nx) + T_cd(2:Nz-1,1:Nx-2)) ...
                        + (1 - 4*fo).*T_cd(2:Nz-1,2:Nx-1);

    i = Nz;
    fo_t = FO(i,2:Nx-1); bi_t = BI(i,2:Nx-1);
    Tp(i,2:Nx-1) = fo_t.*(2*T_cd(i-1,2:Nx-1) + T_cd(i,3:Nx) + T_cd(i,1:Nx-2) + 2*bi_t*T_inf) ...
                 + (1 - 4*fo_t - 2*bi_t.*fo_t).*T_cd(i,2:Nx-1);

    i = 1;
    fo_b = FO(i,2:Nx-1);
    Tp(i,2:Nx-1) = fo_b.*(2*T_cd(i+1,2:Nx-1) + T_cd(i,3:Nx) + T_cd(i,1:Nx-2)) ...
                 + (1 - 4*fo_b).*T_cd(i,2:Nx-1);

    j = 1;
    fo_l = FO(2:Nz-1,j); bi_l = BI(2:Nz-1,j);
    Tp(2:Nz-1,j) = fo_l.*(2*T_cd(2:Nz-1,j+1) + T_cd(3:Nz,j) + T_cd(1:Nz-2,j) + 2*bi_l*T_inf) ...
                 + (1 - 4*fo_l - 2*bi_l.*fo_l).*T_cd(2:Nz-1,j);

    j = Nx;
    fo_r = FO(2:Nz-1,j); bi_r = BI(2:Nz-1,j);
    Tp(2:Nz-1,j) = fo_r.*(2*T_cd(2:Nz-1,j-1) + T_cd(3:Nz,j) + T_cd(1:Nz-2,j) + 2*bi_r*T_inf) ...
                 + (1 - 4*fo_r - 2*bi_r.*fo_r).*T_cd(2:Nz-1,j);

    i=Nz; j=1;  fo_c=FO(i,j); bi_c=BI(i,j);
    Tp(i,j) = 2*fo_c*(T_cd(i-1,j)+T_cd(i,j+1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T_cd(i,j);
    i=Nz; j=Nx; fo_c=FO(i,j); bi_c=BI(i,j);
    Tp(i,j) = 2*fo_c*(T_cd(i-1,j)+T_cd(i,j-1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T_cd(i,j);
    i=1;  j=1;  fo_c=FO(i,j); bi_c=BI(i,j);
    Tp(i,j) = 2*fo_c*(T_cd(i+1,j)+T_cd(i,j+1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T_cd(i,j);
    i=1;  j=Nx; fo_c=FO(i,j); bi_c=BI(i,j);
    Tp(i,j) = 2*fo_c*(T_cd(i+1,j)+T_cd(i,j-1)+2*bi_c*T_inf) + (1-4*fo_c-4*bi_c*fo_c)*T_cd(i,j);

    dT_max  = max(abs(Tp(:) - T_cd(:)));
    T_cd    = Tp;
    time_cd = time_cd + dt;

    T_peak_history(end+1) = max(T_cd(:));
    T_time_history(end+1) = time_cd;

    if time_cd >= nextFrame_cd
        set(hImgCd, 'CData', T_cd);
        set(axC, 'CLim', [T_inf, T_jmax]);
        title(axC, sprintf('Cooldown: t = %.3g s', time_cd));
        drawnow('expose');
        writeVideo(vw_cd, getframe(figAnimCd));
        nextFrame_cd = nextFrame_cd + frame_dt_cd;
    end

    figure(figCool);
    while (snap_k_cd <= num_snaps-1) && (time_cd >= snap_times_cd(snap_k_cd))
        subplot(2,4,snap_k_cd);
        imagesc(x*1e3, z*1e3, T_cd); set(gca,'YDir','normal');
        set(gca, 'CLim', [T_inf, T_jmax]);
        colormap('hot'); colorbar;
        title(sprintf('t_{cd} = %.3g s', time_cd));
        xlabel('x (mm)'); ylabel('z (mm)');
        hold on;
        yline(t_die*1e3,         'c--','LineWidth',1);
        yline((t_die+t_tim)*1e3, 'b--','LineWidth',1);
        hold off;
        snap_k_cd = snap_k_cd + 1;
    end

    pass_ct = pass_ct + (max(T_cd(:)) - T_inf < 1.0);
    if max(T_cd(:)) - T_inf >= 1.0, pass_ct = 0; end
    if pass_ct >= hold_steps
        fprintf('Cooldown converged: t = %.4f s, %d steps, peak T = %.3f K\n', ...
                time_cd, step_cd, max(T_cd(:)));
        break
    end
end

figure(figCool);
subplot(2,4,num_snaps);
imagesc(x*1e3, z*1e3, T_cd); set(gca,'YDir','normal');
set(gca, 'CLim', [T_inf, T_jmax]);
colormap('hot'); colorbar;
title(sprintf('Final (cooled)\nt_{cd} = %.2f s', time_cd));
xlabel('x (mm)'); ylabel('z (mm)');
hold on;
yline(t_die*1e3,         'c--','LineWidth',1);
yline((t_die+t_tim)*1e3, 'b--','LineWidth',1);
hold off;

set(hImgCd, 'CData', T_cd);
title(axC, sprintf('Final (cooled) t = %.3f s', time_cd));
drawnow('expose'); writeVideo(vw_cd, getframe(figAnimCd));
close(vw_cd); close(figAnimCd);
fprintf('Saved cooldown animation: %s\n', videoName_cd);

%% =========================================================================
%  TEMPERATURE HISTORY PLOT
% =========================================================================
figure('Name','Peak Junction Temperature - Cooldown Phase','Color','w');
plot(T_time_history, T_peak_history - 273.15, 'b-', 'LineWidth', 2);
xlabel('Time after cores off (s)'); ylabel('Peak Temperature (deg C)');
title('Cooldown Transient: Peak Junction Temperature vs. Time');
grid on;
yline(T_inf - 273.15, 'r--', 'Ambient (27 deg C)', 'LineWidth', 1.5, ...
      'LabelHorizontalAlignment','right');
yline(T_jmax - 273.15, 'k--', sprintf('Steady-state peak (%.1f deg C)', T_jmax-273.15), ...
      'LineWidth', 1.2, 'LabelHorizontalAlignment','right');

fprintf('\nCooldown complete. Total cooldown time: %.2f s\n', time_cd);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":41.6}
%---
%[output:5e659e70]
%   data: {"dataType":"text","outputData":{"text":"Grid: 11 x 4 nodes | dx = dz = 0.002 m\n","truncated":false}}
%---
%[output:8e3bd423]
%   data: {"dataType":"text","outputData":{"text":"Fo = 0.1091 (lim 0.2425) | dt = 0.003765 s | Bi = 0.03077\n","truncated":false}}
%---
%[output:51d49ccd]
%   data: {"dataType":"text","outputData":{"text":"\nRunning explicit solver (vectorized)...\n","truncated":false}}
%---
%[output:567f3300]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQu4TdX6\/1\/Fli7YhVxK1O8o\/YVyEgoVyaUOIRxOl6NQEhVyjcgluW6JU4dCLilJRUqUcqnoJCeOOEUOTkLZiiOX+D\/fdzWXtcZee68x55pr7jnHfMfzeGrvPeaYc36\/4x3jM8ccc4wCp06dOkWSRAFRQBQQBUQBUUAUCJACBQRgAuSWXKooIAqIAqKAKCAKsAICMFIRRAFRQBQQBUQBUSBwCgjABM4yuWBRQBQQBUQBUUAUEICROiAKiAKigCggCogCgVNAACZwlskFiwKigCggCogCooAAjNQBUUAUEAVEAVFAFAicAgIwgbNMLlgUEAVEAVFAFBAFBGCkDogCooAoIAqIAqJA4BQQgAmcZXLBooAoIAqIAqKAKCAAI3VAFBAFRAFRQBQQBQKngABM4CyTCxYFRAFRQBQQBUQBARipA6KAKCAKiAKigCgQOAUEYAJnmVywKCAKiAKigCggCoQWYA4fPkzYiPvcc89NqRZ89dVXtGHDBrrggguobt26Ccv77bffaM2aNbR9+3a69NJLqU6dOnTGGWfweXft2kWffvop\/3zddddRuXLloteD33\/77bdUpUoVqlq1KhUoUIDWr19PP\/74Y45rrlGjBmVmZvI5Vq9eTWXLlqXatWtTkSJFUrq\/IBwsXgbBJb1rFC\/1dApKLvEzKE4F8zpDCzC33XYb7dixgwAgTtOCBQuoV69e0cMBGfPmzaPChQvHFdm3b1969dVXo7+75557aPDgwbRlyxZq3bo1IciRzjnnHHrrrbeoYsWK9Oyzz9L48eOjx9SsWZNeeeUVuvvuu2nVqlU5Lvm1116j48ePU\/v27aN\/K126NH344Yc5rsfp\/fr1OPHSr87Yvy7x0r5mfj5C\/PSzO8G\/tlACDKDloYceogMHDtDkyZN55MRuwuhN\/fr1uYwXX3yR5s6dSwsXLqQXXniBGjZsGC1u9+7dXD7g5qmnnqI+ffrQ119\/TZ999hlNmDCBjxs1ahQdOXKEnnzySbrvvvs4T6VKlfiYv\/3tb\/TMM89w2W+88QZlZGTQL7\/8wuVnZ2fTAw88wMDz5ptv8rH\/+te\/6N133+URH0BTz549qVatWnG3d+LECVq7di2PHGF05+qrr055JMqufm7lFy\/Fy9i6JHHpVmSlXo7EpjmxmXptSE8JoQSYHj160Ntvv82K4tXPunXrouoCJPAKJlHCKEjRokX5TwAXvLYBnMyYMYNWrlxJGFkBGAEarLRs2TLq3Lkz9e7dmx588EEaO3YsPffcczRz5kwaMWIEw8y\/\/\/1v+vXXX+mqq66iatWqMXh89NFH\/DoJr5y6dOlCH3\/8MS1fvpxhxUqPPfYYgw2A5cILL2QQadCgAf\/DaE6jRo0Sjr60bduW7xl5MPrTtGlTmjRpUnpqWJpLFS\/FS4nLNAeZw+IlNs2JTYdVIO2HhRJgoGpuQ5v\/+c9\/6MYbb0woPEABIyNI27Zt45EWvALCCMnmzZupWbNm1K5dOwYTK+GVUr9+\/Rhc7rjjDgYXjLRg9GXAgAF01llnRQEKsILXPhg9QcJIC4AGCf\/FKyvMg0H6xz\/+QXfeeSfdf\/\/91L9\/\/+j1xF44AAXgg7kxVtq6dSs1btyYR2XGjBlD77\/\/Pr9Gwyut8847L+0VLh0nEC\/FS4nLdERW6mVKbJoTm6nXBvdLEIBR5sBgBAavVxIljLhYk34t0Ln99tspKyuLX8cAUDBHBYBiJbz2wYgMoAZwM23aNBo+fDjPcUE+jLwAIDD0fdlll\/EIC0ZakDCnBdCC11KYnAsYuvbaa\/lvGJUBfHz++ed0\/vnn82TgevXq8YgS5r3MmTOHnn76aRo2bFjcvBi8PkK+PXv2cDk4H0aNWrZs6X7t8qjE3BpJ8dIjA1w8jXhpTlyiWoifZvnpYqi7UpQAjAIweDUEwEiUACJlypThPx06dIjnqAAoABbWq6JHH32UHn744ejh1qslDKfiH0Bm6tSpPPdl5MiR9M9\/\/pNHb\/73v\/\/xKymMjIwePZrhA6+AAEh43YVjrddTFqzcdNNNDERWwisozMvB66CNGzfSn\/70p+gIjZUHnTrgC3NwMIpjvUpbvHgxVa5c2ZVK5XUhuTWS4qXXTqR+PvHSnLjMC2AkNlOPFSmBKNQAgwmveKVzww03ROvCzp07GQISpffee4\/+8Ic\/RP+EERWM1nTt2pWWLFnCoyQYMcGoyhNPPEF33XUXj2wAdPA6B6Mm48aN4\/\/HyAlGVvClEeagHDt2jCEIXyzde++9dMUVV\/BoyiOPPEI4L748ev3113meC2ADkITXPph3YyVrVAajLphDgxGa2FEb5LPABvNk\/vrXv\/LfATEYsVEn+wYlQNDpiZfipVVfJS79E7kSm+a0s\/6pVaevJLQAY32mDJhw+ik1Rk7wyshalwXzUTDfBa+AMHHXGnVRP7fGuTFf5uDBgzxCgpEQJExGBNRgorA16oLfA2Q6dOjAMIOEV1b4h8nDsV9QAb4wQdd6PQQwQr4zzzwzru7h1RXADQn336pVq7jXXn6sqHldk3gpXsbWD4lL\/0SwxKY5semfWiUAwwrgNRDmnqQyeRWL1AEcMA\/F+kIpkdFHjx6l\/\/73v\/xlET6Fjk0ADixkV6pUqbjfY74Kyi5fvnwOCMmtMuF+vvvuO56rU7JkyVzrHCYI47yYA1OwYEE\/1k1b1yReipexFUbi0lb4pDWzxKY5sZnWiuKg8NCOwDjQSg4RBUQBUUAUEAVEAZ8oIACjYQRWt8Wy\/kFJmMuCOS3JUqL70j02Wdl+\/XvQvISOOp6Il36tcfHX5dRL3XoQDBUSX2XQYlPHS9xpGGPTq3ooAKOhNNZn2bF9O11SsSL\/162UrvJQLtapSZZwX9u2\/Ssu26WXXql1bLKy\/fr3oHgJ\/az6oeNnmL2M1cqtepefsZnIS9yXxKYzd\/PTy4hv4WtnnTll\/ygBGA3NrE6veGYmZR84oHGEXpZ0lafT4Z0OrC8UgLkmFACTLu31nNfLZV2jjp+RRjKcXkJNv\/uZqpeReJXY1Iuc+Fzpqhs6cRnWdtaJT06OEYDRUM0CGI2svshiL7A+VgCmXigAxhdGaV6Ejp8RgBEvNSXNt2xOvYx0hBKb+WZcghPreHkaYMIVm175JACjobTZALNUAZhGAjAadcLLLDoNZQRgxEsvfXFyLqdeRjpCiU0nmqfrGB0vTwNMuGIzXZqr5QrAaChtNsC8pQDMnwRgNOqEl1l0GsoIwIiXXvri5FxOvYx0hBKbTjRP1zE6Xp4GmHDFZro0F4BxoKzZAPOqAjBtBGAc1JF0HqLTUEYARrxMpw9ulO3Uy0hHKLHphgdulaHj5WmACVdsuqVxsnJkBCaZQr\/PInfz6yONU6aUxV5gRVaKtNKll94tAJOS+u4frONnBGDES\/fVd7dEp15GOkKJTXfdSK00HS9PA0y4YjM1ZfWPFoDR0MrsEZjJCsB0FYDRqBNeZtFpKCMAI1566YuTczn1MtIRSmw60Txdx+h4eRpgwhWb6dJcLVcARkNpswFmnAIwjwnAaNQJL7PoNJQRgBEvvfTFybmcehnpCCU2nWiermN0vDwNMOGKzXRpLgDjQFmzAWaEAjD9BWAc1JF0HqLTUEYARrxMpw9ulO3Uy0hHKLHphgdulaHj5WmACVdsuqVxsnJkBIaIfvrpJ1q6dCnt3r2bWrZsyRscxiazAWaQcq9DAw0wpnkJc3QaygjAmOUl7j0vP4MWl6l4GekIJTaTdWhe\/l0nLk8DjHmx6aXWuZ0r9ACDHZ9bt25NxYsXp2rVqhG2f\/\/ss8\/idnIOWkNpL7AeUwBmXA6AOXbsGI0bN47Wr19PNWrUoE6dOlFmZiatWbOGpk+fzjtfd+\/enSpUqBAta+DAgXTfffflgMF0VnoTvbTX6SX3EuWZ4mfQ4jIVLyMdocRmOtsPu2WHtZ21q1M684ceYDDyMnr0aHr\/\/fdZ53feeYeuuOIK3r\/CSkFrKO0FVlcFYCbnAJjZs2fTihUrqHfv3jRz5kwqXbo0tW3blpo2bUpDhw6l7du38wjWwoULaePGjZSVlUXLly+nxYsXU+XKldNZf+PKNtFLe51eci9Rnil+Bi0uU\/EyAjASm541JhonCms7qyGNZ1lCDzBTpkyhDz74gAoVKkTHjx\/n0Rh0zrEpaA2lvcDqqNzrizkAZuXKlVSuXDkeTRk1ahQVKFCAR2Lmzp1L06ZNo1OnTlHVqlVp7dq1dPjwYdqxYwc9\/vjjNGnSJE8BxkQv7XV6yb1Eeab4GbS4TMXLCMBIbHrWM2qcKKztrIY0nmUJPcAMGDCAO2KMLGBovU+fPvyEWqlSpbgRmIPZ2a5u5Jguh7FxWbHixbXmsUTmTXSIXkpW1leUlfXPXI+tU6cO7dmzh2bMmMF5MGeof\/\/+fDxGY8aPH0+XX345\/9yiRQsaOXKkpwBjmpfQUddPu16i7KD7aY2SmhabqpfwSmIzXS2ms3J14zICnma1s84US89RoQcYjCj88MMPPMcDqW\/fvnTJJZfQgw8+GAcw5i5k10oZgXk9B8Ds2rWL5wRlZGTQokWLaMKECdSrVy9avXo1DRs2jI9v0KABLViwgIoVK5ZvAGOil\/ae2pN7ifJM8dPsEZh4LyMdocRmerpBZ6XaG4FJHptBiUtnaqXnqNADzJIlS3ikAB3zmWeeySMHzzzzDF199dUhAZhbFYB5LwfAAOpq1arF2qxbt47Gjh3L84Y6duxI0G\/nzp3UoUMHntRrpfwYgTHRS3sAk9xLC9JN8NNsgIn3MgIwEpvp6QadlWoPYJLHZlDaWWdqpeeo0APMb7\/9RmPGjKH58+fzPJjmzZvz\/A3M87BS0BpKe4FVXwGYj3IADCbmIriQTp48SfjCqHbt2tSvXz8Gmr179zLUNGrUKF8BxkQv7QFMci9Rnil+Bi0uU\/EyAjASm+npBp2VGtZ21pla6Tkq9ABjyXro0CEegSlSpEgOpYPWUNoLrOsUgPks1zkw+\/bti\/u8HAfid\/iMOpFu6amyyUs1yUt7nZ6+l5Z3eDUYm4LkZ9DiMhUvIwAjsZk8+r3LEfZ21julcz+TAIyGC0FrKO0FVjUFYDZoTQDWkM2XWYLmpb1OT7z0ZaWLuSid2IxM+oz3MgIwEpt+8lfHy4hvOf003UuvfBKA0VA6aJ2evcC6TAGYbwVgNOqEl1l0\/Iw0kuKll744OZdTLyMdocSmE83TdYyOl6cBJlyxmS7N1XIFYDSUNhtgyikAs1sARqNOeJlFp6GMAIx46aUvTs7l1MtIRyix6UTzdB2j4+VpgAlXbKZLcwEYB8qaDTAXKADzowCMgzqSzkN0GsoIwIiX6fTBjbKdehnpCCU23fDArTJ0vDwNMOGKTbc0TlaOjMAkU+j3d5jmrgNzjgIwhwVgNOqEl1l0GsoIwIiXXvri5FxOvYx0hBKbTjRP1zE6Xp4GmHDFZro0lxEYB8qaPQJzpgIwvwnAOKgj6TxEp6GMAIx4mU4f3CjbqZeRjlBi0w0P3CpDx8vTABOu2HRL42TlyAhMMoWMH4E5ogBMEQEYjTrhZRadhjICMOKll744OZdTLyMdocSmE83TdYyOl6cBJlyxmS7NZQTGgbJmj8AcUAAmUwDGQR1J5yE6DWUEYMTLdPrgRtlOvYx0hBKbbnjgVhk6Xp4GmHDFplsaJytHRmCSKWT8CMz3CsCUEYDRqBNeZtFpKCMAI1566YuTczn1MtIRSmw60Txdx+h4eRpgwhWb6dJcRmAcKGv2CMx2BWAqCsA4qCPpPESnoYwAjHiZTh\/cKNupl5GOUGLTDQ\/cKkPHy9MAE67YdEvjZOXICEwyhUwfgfmXElhXSiOpUSU8zaLTUDLAiJee+uLkZE695I5QYtOJ5Gk7RsfLKMCELDbTJrpSsACMhtJGj8BsUACmmgCMRpXwNItOQ8kAI1566ouTkzn1kjtCiU0nkqftGB0vowATsthMm+gCMPalNRpg1ikAc21OgDl27BiNGzeO1q9fTzVq1KBOnTpRZmYmrVmzhqZPn86bOXbv3p0qVKhAR44coYkTJ9KmTZuoYcOGdPfdd8cJ\/tNPP9FLL71E\/\/jHP\/jv9957L51xxhn2TXF4RNC8xG3qNJQMMBpeojxT\/AyTl9wRSmw6jPr0HKYTl1GA0YhNU+IyPWonLlVGYDTUDlpDaSuw1igAUycnwMyePZtWrFhBvXv3ppkzZ1Lp0qWpbdu21LRpUxo6dCht376dli5dSgsXLqSxY8fS1q1bqWPHjjRw4EAaNGgQ1a1bN6rymDFj6LfffqM2bdrQ448\/Tg8++CDdfPPNGi64kyVoXtoCGA0vUZ4pfobJS+4IJTbdaQRcKkXaWZeETKEYARgN8YLWUNoKrI8UgKmfE2BWrlxJ5cqVo4oVK9KoUaOoQIECPBIzd+5cmjZtGp06dYqqVq1Ka9eupVatWtGIESOoevXqNGPGDMrOzqYePXpEVX7yySfpvPPOo\/bt29MjjzxCXbp0EYBJUgd1\/OQRGA0vcSpT\/AxaXNqCUcVLBhiJTY3W2rssOnEZHYHRiE1T4tI7B4gEYDTUDlpDaSuw3jsNMFmzMinr5eK5foVUp04d2rNnD4PJtm3baPfu3dS\/f39WEKMx48ePpyZNmvCrpmLFitGyZct4VGbSpElRlT\/88EO677776JxzIktrv\/\/++zyi41UKmpe2Oj0bXqLcoPsZFi\/hlcSmVy2E\/nmkndXXKl05BWA0lA1aQ2krsBYpIzC35RyB2bVrF5UsWZIyMjJo0aJFNGHCBOrVqxetXr2ahg0bxgo2aNCAFixYQC1btqRZs2ZRmTJlaMmSJbRhwwbq27dvVGXk69evH4+6jBw5ks466yzq2bOnhgvuZAmal7YARsNLlGeKn2Hykp\/kJTbdaQRcKkXaWZeETKEYARgN8YLWUNoKrAUKwLTMCTAAkFq1alGLFi1o3bp1PM9l9OjRPM8FkLJz507q0KEDT+pF3nr16lHjxo15dKZmzZrUrFkzOnjwIJUqVYpuu+02Gj58OFWrVo1fP+EVkwBM3pVQx09+haThJc5kip9Bi0tbMKp4yQAjsanRWnuXRScu2TfN2DQlLr1zQF4haWkdtIbSVmDNUwCmbU6A2bhxY3QU5eTJkzw5t3bt2jySAqDZu3cvQ02jRo1o8+bN\/OVR0aJFedQGr5vwO0DK8uXL6YMPPiDMg7nwwgtZe7x2uuiii7R8cCNT0Ly01elpeInyTPEzTF5yRyix6UYT4FoZ0s66JqXjgmQERkO6oDWUtgJrhgIw9+S+Dsy+ffsYSmITfofPqIsUKRL99fHjx2n\/\/v38GilRAgTh7xiR8ToFzUtbAGPDS5QbdD\/D5CUDjMSm181FnueTdjb\/7TAeYPDq4osvvuBRAkwsrVSpEpUvX54KFiyYQ328DilcuHCOr2KC1lDaCqypCsDc79+F7MLopS2ACZCXuK9U\/QxaXKbiJQOMxGb+95gxV2BqO+srkZNcjLEA8+233\/JrjXfffZcluOyyywhP\/lizBAnrkGDxtbJly\/LPW7Zs4S9o8Dt83hubgtZQ2gqsyQrAdPUfwITZS1udXgC8xP245WfQ4jIVLxlgJDZ91bea1s76SlzNizESYObPn8+rwWKNESyiBkixRlxOnDjBnwJjtOXll1+mjz\/+mFePxfolmLeByaqhApgsBWB6+Atgwu6lrU7P517iXtz002iAUbxkgJHY1OzWvMlmC2ACEJveqObuWYwEGIyyXHzxxQlfE8XKd\/ToUX5lhEmlJUqUIMANUiKAOZidTdkHDrirfhpKK56ZScWK576WizqytG2MAjC9\/AUwYfYSXun6yV86+NxL3I+bfuKekUyLzUReMsBIbKahxXRWpG5csm8BiU1nSuTvUUYCTKyk+LQXy+D\/\/PPPcUoPHjyYJ55isbXnn3+eV5W1FlwL1QjMMAVgBvoLYMLupa0RmAB5iftKNTaNHoFRvOSOUGIzf3tL5ey2RmACFpu+EjqPizEaYPC0h4XTbr\/9dl4KPzZheXuMvmBPH3wRg896v\/vuO86ClWLxz0pBayhtBdYgBWCG+hNgwuqlLYAJiJfWSEyqsRm0uEzFSwYYiU1f9asmtrO+EljjYowGmFWrVtGQIUN4ufrcElYlxaskJGxUiH19unbtGre8fdAaSluB1U8BmJH+BJiwemmr0wuIl7gnN\/wMWlym4iUDjMSmRpfmXRYT21nv1HPnTEYDDD7TbNiwIT3xxBNUv359XrbeShh9URNeIWEeTKheIfVSAGaMPwEmrF7a6vQC4iXuyQ0\/jQYYxUsGGIlNd3o9l0qxBTABik2X5PGkGKMBBjBy55138n48asLvsCuyTgpaQ2krsLopADPJnwATVi9tAUxAvMQ9ueFn0OIyFS8ZYCQ2dZprz\/KY2M56Jp5LJzIaYLArMj6PnjdvXo4dj7F8fYECBbRkDFpDaSuwOisA84I\/ASasXtrq9ALiJe7JDT+DFpepeMkAI7Gp1V57lcnEdtYr7dw6j9EAg4mfWLDuk08+SfpJdV6CXnjhI7R\/7z5NzeckyNde81hkS\/14W4H1VwVgXvInwLjlZaTTq63pR+peJD6Rbn2InF\/HT\/5UMyBe4p7c8NOel+7EVk4\/db2MnN+plwwwhsemtLOazZJkiypgNMAcOHCAmjdvzgvUYXfk2H187rjjDsrIyNCqCgisfXsjmw8mSwWof44sp2hEssNOG+LC8TqNJDeI6PT+ogDMLH8CjFte4p6\/295Jy490eIkT69YH6\/w6fgbJS2jghp92vMQ50+GnrpfW+Z16yfFqeGxKO6vVLEmmGAWMBpidO3fyBN5EacqUKXEbEOZVK4I2VK3TSEYBpq0CMPNyAsyxY8do3LhxPOxfo0YN6tSpE2VmZvI6HtOnT+fNHLEFQ4UKFXhVY6yCvGnTJp5AjZ2pYxO+8vrb3\/5GS5cu5VWP8cWXzlyksHoJ7XT8ZIDR8BLlmeJn0OIyFS85XiU2fdV568Rl0NpZXwmscTFGA4zG\/WtlCVpDaSuwmisA82ZOgJk9ezYvBti7d2\/+1Lx06dK8fk7Tpk1p6NCh\/DoAQLJw4ULef2rr1q3UsWNHGjhwIA0aNIi3c7AS5iNhGwf8\/oUXXqArr7wyB+RomeIwU9C8tNXpaXiJ8kzxM0xeckcosekw6tNzmLSz6dHVTqlGAww2b3zuuee4w8zOzo7TBR0uRg50UtAaSluB1UQBmCU5AWblypW8EGDFihVp1KhRPPkZIzFYvXjatGm8dk7VqlVp7dq1PGl6xIgRVL16dZoxYwbrjkUDrXTvvffSbbfdRueccw5Vq1aNLrjgAl5QMFkKq5e2AEbDS5Rnip9Bi8tUvGSAkdhM1kx4+ncT21lPBXThZEYDzFdffcVzYNChlilTJk6u66+\/Xntib9AaSluB1fA0wGRty6Ssbbnvo1SnTh3eCBNgsm3bNtq9ezf17x+Z84PRmPHjx\/OO3njVVKxYMd6mAaMy1hYNyIcysIVDvXr1aNGiRfy36667LmlVDquXtjo9G15aXgTZz6DFpVMvcZzEZtImwvMMJraznouY4gmNBph169ZRz549ecfpVFLQGkpbgVVfGYH5KOcIDFYrxgRoTHoGdEyYMIF69epFq1evpmHDhrG0WBZ+wYIF1LJlS5o1axYDI0a+sN5O3759o\/IDYEaOHMkLC77zzjv0wQcf0JgxY5LaE1YvbXV6Gl6iPFP8DFpcpuIlj8BIbCZtJ7zMYGI766V+bpzLaIDBFgHt2rXjeRboYGNX4q1Zs6b2CMyFjzxC+\/dpfkY9J8Gnt+1tfGrpwvG2Aus6BWA+ywkwABBMuG3RogUBJDDPZfTo0TzPBZCCCbYdOnTgSb3Ii9EVfPWF0Rno3KxZM155tVSpUjwBGKM1+ArsxRdfpMOHD9PDDz+ctC675SV3erU1P6N2wYuEN6ZbH34\/v46fPIlXw0tcjyl+2vISN54OP3W9\/P38Tr1kgDE8NqWdzf92NmlD7LMMRgMMOscbb7yRfvzxxxwL2dmZA5M5ZQodvPhiLesK3H57jnyn3n5b61hkcuN4nUaSG0R0elcrALM+J8Bs3LgxOoqCuSiYnFu7dm3q168fA83evXsZaho1akSbN2\/mSbn4dB2jNnjdhN9hJGz58uW0ZcsWPh5fK5UoUYJf75UtWzapPm55yZ\/eTpyY9HxueZHoRLr1waoLOn7qeonrMcVPO16my09dL63zO\/WS49Xw2JR2Nv\/bWa2G0UeZjAaYTz\/9lJ\/usXGczkTR3HwJ2lC1TiMZBZgqCsBszH0dmH379sWtpYMy8DtMhsa8FisdP36cd\/hW5x3F6ot1QPAptm4Kq5fQR8dPBhgbXlrexa6NFDQ\/gxaXqXjJ8SqxqdtceJJPJy6D1s56IpyLJzEaYH766Seea4F5FmpDbUfDoDWUtgLrcgVgtvhzIbuwemmr0wuIl7gnN\/wMWlym4iV3hBKbdprttOc1sZ1Nu2gun8BogMFXMtjMEV9aYC5G7KuK4cOHy0J2eGqvqADMdn8CTFi9tNXpBcRL3JMbfhoNMIqXDDASmy53f6kVZwtgAhSbqani7dFGAwxeY+DLmETpnnvu0X6tFLSG0lZglVEA5nt\/AkxYvbQFMAHxEvfkhp9Bi8tUvGSAkdj0tndMcjYT21lfCaxxMUYCzIcffkiXXXYZlS9fPk8JsF7J1VdfnVSmoDWUtgKrhAIw+\/0FMGH30lan53MvcS9u+hm0uEzFSwYYic2kbbWXGUxqZ73Uzc1zGQkwn332GS97\/8c\/\/pFuuukmuuiii\/gf9oDZsWMHL33\/1ltv8c9vvPFGUj2D1lDaCqxiCsAc9BfAhN1LW52ez73EvbjpZ9DiMhUvGWAkNpO21V5mMKmd9VI3N89lJMBAIHym+9prr\/FqsBhpwWe4VsLePK1bt+b1Sc4+UxBjAAAgAElEQVQ444ykegatobQVWGcrAPM\/fwFM2L201ekFwEs3\/QxaXKbiJQOMxGbSttrLDKa1s15q59a5jAUYVSCsVQKIwUReu59UB62htBVYZygAc9J\/ABNmL211egH0EvfnNDaDFpepeMkAI7HpVr\/nSjmmt7OuiJTmQkIDMKnoGLSG0lZgnVAApqD\/ASZMXtrq9MTLVKqGJ8fqxCav6aN4yQAjsemJR7on0fGSfUvgp+le6mqYaj4BGA0FjQaYIwrAFBGA0agSnmbRaSi5kRQvPfXFycmceskdocSmE8nTdoyOl1GACVlspk10pWABGCLaunUrL3N\/3nnn8R4+WOI+NhkNMD8rAFM02ABjmpe2RmAM8xL3npefQYvLVLzkjlBi06t+Ues8tgDGwNjUEinNmYwGGEzkxUaDlSpVipPxo48+ouuvv543c8TfsbkgVuzF\/JipU6fyLsuxy+AHraG0FVj7FYAp4U+ACauXtjq9gHiJe3LDz6DFZSpeMsBIbKa5O7RXvIntrD0F8j+30QCDTQTxpdFdd91FTzzxRHT3aTR8GzZs4BGXefPm8S7KWVlZ7AY2ImzYsCH\/10pBayhtBdYuBWAu8ifAhNVLW51eQLzEPbnhZ9DiMhUvGWAkNvO\/x4y5AhPbWV8JrHExxgNMmzZteNPAChUq0LPPPkvFihXjSVUWwGBF0FOnTvFeSfh\/rBszc+bMuAXukP9gdjZlHzigIWn+ZimemUnFihenbdu2Jb0QnjexXQGYijkBBuvljBs3jj9Hr1GjBnXq1Ik1BfhNnz6dN3Ps3r07a4wn64kTJ9KmTZtygGDsBaGsadOm0aRJk5Jep9Xhhc1L3Leun7peokxT\/MQ9I5kWm4m8ZICR2NRqK7zIpBuX7FuA2lkvtHPzHMYDTLdu3Wj+\/Pn00EMP0a5du+ill16iW265JQowlphvv\/02DR06lFq1akV9+\/aN0zhoT3p2ngy2KABzeYJGcvbs2bRixQpeHBBwV7p0aWrbti2\/eoNmWBhw6dKltHDhQho7dizPW+jYsSMNHDiQBg0aRFh3Jzb98ssvPDIG8HnnnXe06jOe2MPopZ2ndh0vUZ4pfgYtLlPxEsdKbGo1FZ5lMrGd9Uw8l04UCoDBBN3jx48TNnB8\/fXXeT0YawTmxIkTNHjwYP552LBhVL169RzSBq2htBNYXykAc1UCgFm5ciWVK1eOKlasSKNGjaICBQrwSMzcuXN5FAUjWFWrVqW1a9cyAI4YMYJ1nDFjBmVnZ1OPHj3iNH3ssceoWrVq\/PrOLsCEzUs7nZ6OlyjPFD+DFpepeIljJTZd6vVcKsbEdtYlaTwrxmiAwY636ET79+8fFRRPn08\/\/TRP1C1atChv9vjiiy8y2BQqVIjzoYPGPysFraG0E1jrYgDmhcxM+nser5\/q1KnDO3tDU7yigr6WthiNGT9+PDVp0oRfNeFVHVZBxqhM7Gsi6P3FF1\/Qgw8+yK+idAEmrF7a6fTseIlyg+5n0OLSqZc4TmLTsz5R+0QmtrPaN++TjEYDjI7G6IBfeeWVuKy9evWirl27hgJg1igjMHUSjMDg1RvmCGVkZNCiRYtowoQJBI0AgRi1QmrQoAHDYMuWLWnWrFn8FdeSJUt4ZMt6JXf06FGqXLkyz5fByBdGYKD\/HXfcoWNV0jwmemmn09PxEuWZ4qfJAKN6ycApsZm0DfAygx2A0YnNoMSllxonO1foASaZQPh70BpKO4H1kQIw9RM0kgCQWrVqUYsWLWjdunU8z2X06NE8zwWQgk\/RO3TowJN6kbdevXq8ng6AombNmjzf5eDBg3T++efTe++9x5Lj1dJzzz1HzzzzDN1www06NriSJ2he2gEYHS9Rnil+hslL+Cax6UoT4Foh0s66JqXjggRgNKQLWkNpJ7DeVwDmlgQAs3HjxugoysmTJ3lybu3atalfv34MNNjLBlDTqFEj\/jwWn6Dj9RxGbfC6Cb\/r2bMnLxZoJbyKAgDpvkLSsEkrS9C8tAMwOl6iPFP8DJOX8E1iUyvEPcsk7axnUud6IgEYDQ+C1lDaCax3FIBpmgBgLIn27dvHUBKb8Dt8TVSkSJHorzFhGp+kxy4GqCGzJ1mC5qUdgLHjJcoNup9h8hJ+SWx60kRon0TaWW2p0pZRAEZD2qA1lHYCa4ECMC3zABgNqXyfJWhe2gEY8dL31Y90YhN1VPUSdyax6S9\/dbzEFSfy03QvvXJKAEZD6aB1enYCa54CMG0FYDRqhLdZdPxEHRUvvfXFydmceolzSWw6UTx9x+h4aQFM2GIzfarHlywAo6G0yQAzSwGYvwjAaNQIb7PoNJSoo+Klt744OZtTL3EuiU0niqfvGB0vLYAJW2ymT3UBGNvamgwwLyoA01EAxnb9SPcBOg0l6qh4mW4nUi\/fqZc4s8Rm6vq7WYKOlxbAhC023dQ5r7JkBEZDaZMBZrICMF0FYDRqhLdZdBpK1FHx0ltfnJzNqZc4l8SmE8XTd4yOlxbAhC0206e6jMDY1tZkgMlSAKaHAIzt+pHuA3QaStRR8TLdTqRevlMvcWaJzdT1d7MEHS8tgAlbbLqps4zApKimyQAzRgGYXgIwKdYW9w\/XaShRR8VL97V3u0SnXuI6JDbddiO18nS8tAAmbLGZmrL6R8srJA2tkgHMb6dO5SjlzJi9lBKdQj3Gzfx2AmuEAjD9BWAond6gLtgtX8dP1FHxMmek2dU61fy4grxi2amXKDfssSntrEZnFbIsAjAahpsMMEMUgBksAGMbMFLt9JLBq26nJ16aAzCql7izsMdmkAEmbLGp0a26kkUARkNGkwGmnwIwIwVgAgsw4qU5AKN6iTsLe2wGGWDCFpsa3aorWQRgNGQ0GWB6KQAzRgAmsAAjXpoDMKqXuLOwx2aQASZssanRrbqSRQBGQ0aTAaa7AjATBWACCzDipTkAo3qJOwt7bAYZYMIWmxrdqitZBGA0ZDQZYB5UAGZKAoA5duwYjRs3jtavX081atSgTp06UWZmJq1Zs4amT5\/Omzl2796dKlSoQEeOHKGJEyfSpk2bqGHDhrwzdWz65Zdf6Omnn6ZvvvmGd6++8847eedqr1IyL51MsvXLHBgdL3F\/pvjpRy\/dmsSreolywx6bQQYYndg0JS69astxHgEYDbWTNZRBDqy\/KgDzUgKAmT17Nq1YsYJ69+5NM2fOpNKlS1Pbtm2padOmNHToUNq+fTstXbqUFi5cSGPHjqWtW7dSx44daeDAgTRo0CCqW7duVOWnnnqKDh48SF26dKGsrCyqWrUqde7cWcMFd7Ik8zLIAKPjJe7PFD\/96KVbAKN6iXLDHpvSzgannXWntU5eigBMco14N9EdSkcfe1iQA+svyn3NSgAwK1eupHLlylHFihVp1KhRVKBAAR6JmTt3Lk2bNo1OnTrFILJ27Vpq1aoVjRgxgqpXr04zZsyg7Oxs6tGjR1QuQEvr1q25PMDQli1baPjw4RouuJMlmZdBBhgdL3F\/pvjpRy\/dAhjVS5Qb9tiUdjY47aw7rXXyUgRgkmtkNMC0jQGYTZmZtLF4cdq2bVtCVerUqUN79uxhMEGe3bt3U\/\/+\/TkvRmPGjx9PTZo04VdNxYoVo2XLlvGozKRJk3KUhxGbAQMG0EsvvURVqlTRcMGdLH7s9Nz6jNqOl1Az6H760Uu3ACbWS5QpsZlz\/aRkWqf7YURneQNcA+qpndgMely601LrlSIAo6FTsoYyyE8GdygjMG8kGIHZtWsXlSxZkjIyMmjRokU0YcIE6tWrF61evZqGDRvGCjZo0IAWLFhALVu2pFmzZlGZMmVoyZIltGHDBurbt29UZcyRwauoQ4cO0eDBg3lUx8uUzMt0N3pOytdpKHFfOl7i\/Kb46Ucvk3WqTr1EuWGPTWlng9POetWmC8AQ8UgCRgsKFixIt956K5UoUSJO\/2QNZZAD6zYFYBYlABgASK1atahFixa0bt06nucyevRonucCSNm5cyd16NCBJ\/Uib7169ahx48Y8OlOzZk1q1qwZz3spVaoU\/x1zaB555JG01PFUvXQCGH6ZxKvjJe7PFD+TxWV+eOkWwKheotywx6a0s\/5pZ9PSeDsoNPQAg44VX8tcc801DC6Yx7F48WIebbBSsoYyyIHVUAGYZQkAZuPGjdFRlJMnT\/Lk3Nq1a1O\/fv0YaPbu3ctQg6+KNm\/ezF8e4csijNrgdRN+17NnT1q+fDlde+219Ouvv9JZZ50VHbnBV0luJDe8zI9Oz61XSDpe4v5M8TNZXOaHl24BjOolyg17bEo764921o222q0yQg8wmEiK0Rf8FwmjB5MnT46bl5GsoQxyYNVXAOajPNaB2bdvH0NJbMLv8Bl1kSJFor8+fvw47d+\/n18jeZnc8DI\/Oj23AMaOl7jPoPuZLC7zw0u3AEb1EuWGPTalnfVHO+tlm57sXKEHGMzDOPvss6lPnz6sFT7pxUgCvpQJwwhMbQVgPgnwQnZueJkfnZ5bAGOSl\/AhmZ8mA4zqJfQIe2wGGWBMi81kYOHV30MPMN26deN5GtaCa\/gy5pJLLolbmwQN5cHsbMo+cCChL34KrOKZmVQsjy+JYm8A9\/VHBWA+DzDAuOGl3wBG10\/TvIQPyfzEPSPZic1ksJjqfKa8RmBS8RLlhj02pZ31CguCc57QA8yUKVPo559\/jo7AYIG1du3a0U033RSKEZgqCsBsDDDAuOGl3wAG16P75YpJXuK+k\/lp8giM6iX0CHts+glgdOMS+VBPTYtNvyBO6AFm1apVNGTIEJ64i89LmzdvzqvKxs7fSNZQBjmwLlcAZkuAAcYNL4MMMCZ5CR+S+ZksLvPDy7xGYHQ7PdyX6iWODXtsSjvrF2zwz3WEHmAw4RRD1fgE+PDhw\/xlDfb6UV+1mLoSb0UFYLYHGGDc8DI\/Or1krzV0R2BM8hI+JPPTZIBRvYQeYY\/NIAOMabHpF4QJPcBYRuBTYEzmxRc1akrWUAY5sC5SAGZXgAHGDS+DDDAmegk\/covNZHGZH166NQKjeolywx6b0s76BRv8cx0CMBpeJGsogxxYFyoA84MBAJOXpcm8zI9Oz60RGPEyp\/OpTspN5o3d2NcdTVO9xJ2FPTbtap3uWNbxEteANidssanRrbqSRQBGQ8ZknV6QA6uYAjAHBWDI604vWSep01CijoqX5gCM6iXuLOyxKe2sRmcVsiwCMBqGmwwwZysA8z8BmMACjHhpDsCoXuLOwh6bQQaYsMWmRrfqShYBGA0ZTQaYggrAnBCACSzAiJfmAIzqJe4s7LEZZIAJW2xqdKuuZBGA0ZDRZIA5pQBMAQGYwAKMeGkOwKhe4s7CHptBBpiwxaZGt+pKFgEYDRlNBpgjCsAUEYAJLMCIl+YAjOol7izssRlkgAlbbGp0q65kEYDRkNFkgPlZAZiiCQDm2LFjNG7cOFq\/fj3VqFGD18nJzMzktXOmT5\/On553796dKlSoQEeOHKGJEyfSpk2beJdva4sGDZk9yZLMy3R\/ueCkfN1JvDpe4vym+OlHL6FvXpOynXqJcsMem0EGGJ3YNCUuPWnIfz+JAIyG2skayiAH1n4FYEokAJjZs2fTihUrqHfv3rxrd+nSpalt27bUtGlTGjp0KG3fvp1XL164cCGNHTuWtm7dSh07dqSBAwfSoEGDqG7duhoqe5MlmZdOACPdXy3pdno6XuL+TPHTj166BTCqlyg37LEp7Wxw2llvWnMiARgNpa2GEpuxJdrQ0WlgZWVlUY8ePfgKkn1Kq9NJWten0+HhnLiv7xWAKZMAYFauXEnlypWjihUr0qhRo6hAgQI8EjN37lyaNm0anTp1iqpWrUpr166lVq1a0YgRI6h69eo0Y8YMys7Ojt6jhtRpz5LMS6cA47aXuA47fup6iXJN8TMWYHRjUzfOLD9188dW3ETHpOolyg97bEo7G5x2Nu0NuYzA6Evcvn17+vTTT\/UPyOectWrVojlz5iS9ikT3ldexderUoT179jCYbNu2jXbv3k39+\/fn82A0Zvz48dSkSRN+1VSsWDFatmwZj8pMmjQp6bV4lSFoXkIXHT\/teolyg+5nmLxMVg+C7iXuL2h+6sRlbvdlejvrVXsuIzBeKR3g82CTy5IlS1JGRgYtWrSIJkyYQL169aLVq1fTsGHD+M4aNGhACxYsoJYtW9KsWbN4M8wlS5bQhg0bqG\/fvgG+e\/MuXfw0x1PxUrwMczsrAGNO\/U\/bnQBA8MTQokULWrduHc9zGT16NM9zQfDs3LmTOnTowJN6kbdevXrUuHFjHp2pWbMmQ40k\/yggfvrHi1SvRLxMVUH\/HC9e2vdCAMa+ZqE7YuPGjdFRlJMnT\/Lk3Nq1a\/PO3QAabLYHqGnUqBFt3ryZvzwqWrQoj9rgdVPhwoVDp5mfb1j89LM79q5NvLSnl59zi5f23RGAsa9ZaI\/Yt28fQ0lswu\/wGXWRIkWivz5+\/Djt37+fXyNJ8q8C4qd\/vbF7ZeKlXcX8m1+81PdGAEZfK8kpCogCooAoIAqIAj5RQABGwwh8bYMvagoWLEi33norlShRQuOoxFkwMvHOO+\/ELfDmtHyst7J8+XI677zzeM6JdV1Oy3N8UwE70C19xMv8N96PXkIViU37dcMtL3FmiU37+gfxCAGYJK4dPHiQV5S95pprGBCw1snixYv5ixw7CWulfPnll7xuCj5BBsQgOS0fE2fx6XL9+vWpbNmyNHXqVP4q6Oyzz3bleu3cW5DyOtU79h7FS3847kcvoYzEpv364YaXOKvEpn3tg3yEAEwS97DyLEZf8F8kfGEzefJkqlKlii3fsUx0t27dCO83jx49GgUYp+XPmzePv\/rBgltImDgL0EJy43pt3VyAMjvVO\/YWxUt\/GO5HL6GMxKb9+uGGlzirxKZ97YN8hABMEvcGDx7Moxp9+vThnJ07d+avbVq3bu3Id0AH1k6xRmCclo8hUjxtYFIt\/v+mm25iyMLCcW5er6Ob9PFBTvVOdEviZf4a7UcvoYjEpv164aaXOLvEpn0PgniEAEwS1zBqgrVMrE0JBwwYQJdccgmDjJOkBlaq5b\/99tu8HxGW8Mc6AqmW5+SegnSMm\/qIl\/nrvJ+9hDISm\/r1w00vEwFMquWLl\/peeplTACaJ2lOmTKGff\/45OgLTpUsXateuHY94OElqp+e0\/BMnThCeWrDSLUZ0sPcQktPynNxLEI9xUx\/xMn9rgB+9hCISm\/brhZteJgIYp+WLl\/a99PIIAZgkaq9atYqGDBnCE3exbHfz5s1552Wna5yonZ7T8rFs\/4svvkivv\/46FSpUiO8CmyxiIq+b1+tlZfTiXE71TnRt4qUXjuV+Dj96iauV2LRfL9z0MhHAOC1fvLTvpZdHCMAkURuLsmH4EZ3V4cOHefXZTp06OfZI7fSclo9l+l955ZW468D+RLg2N6\/X8Y369ECneusAjNOyxUtnlcWp3un0EmWLn\/b9dNPLRADjtHzx0r6XXh4hAKOpNpbLx+RYrDqbjuR2+W6Xl457zs8y06mP22W7XV5+6p6Oc6dTn3SUnY4y06FrfpSZbm3cLt\/t8vJD8yCfUwAmyO7JtYsCooAoIAqIAiFVQAAmpMbLbYsCooAoIAqIAkFWQAAmyO7JtYsCooAoIAqIAiFVQADGA+OXLFnC+xXdcMMNrp0NW69v2rSJ2rZt61qZUlByBcTL5BoFJYd4GRSn9K5T\/NTTyaRcAjBpdvPAgQPUrFkzevfdd6lo0aKunQ3bEWAtmtmzZ1PFihVdK1cKyl0B8dKc2iFemuMl7kT8NMtP3bsRgNFVKiYfPslbv349VahQgUqVKsULV1k\/Y2n\/2DR27Fje+wif42GTN+TFFgDffPMNXXvttXTmmWfyBpH4uum6667jtVwwunLhhRfyujOHDh3ilYCxPDkWrQOsVK5cmU\/x97\/\/nTZv3kzjxo1zcBdyCBQQL82pB+KlOV5KbJrlZbruRgDGobJYkfeHH36g+fPn84JygAkscJeZmRlX4i233EKDBg2iunXr0jPPPEOvvfYa\/71IkSJ05MgRHpU5\/\/zz6euvv6aHHnqIHnjgAWrRogVha\/nixYvTnj17eOuCHTt20JVXXknr1q3jMmrUqMHwgtEdvEpCeZKcKSBeOtPNj0eJl350xfk1iZ\/OtQvDkQIwDl3GrtI333wz\/elPf6K5c+fyRorqHBeMvGC0ZPny5TxyAoB5+eWX6fPPP+cn\/6pVq9Kjjz5KDz\/8ME2bNo3ef\/99XpwOAIORHEDRypUr6Z577qE5c+ZQrVq1eKG6atWq8WJ1GJ1BGdj9FqM5kpwpIF46082PR4mXfnTF+TWJn861C8ORAjApuIwNvnr06EF\/\/vOfafjw4TlK2rp1KzVu3JhHVzIyMhhg8Brp2Wef5byXXnopbwVw9dVX88Zvzz\/\/PC1atIgBpk2bNtS+fXt+nQRIQlkFCxYkbCaJCcHYuBEJ4ILNHJs0aZLCncih4qU5dUC8NMdL3In4aZafbt6NAEwKagImMPpy2WWXMXgULlw4rjQLPr788kt+VQSAwWsjbMJoAcxbb71FVapUyQEwGDoFlFhlbNu2jY9RAeaqq66iqVOn8vwZSc4VEC+da+e3I8VLvzmS2vWIn6npZ\/LRAjAO3f3www\/pvvvuozfeeIMefPBB3uSxT58+caVhF2vsEv3mm28SQMNtgDl48CCP3vzjH\/\/IMffG4W2F8jDx0hzbxUtzvMSdiJ9m+en23QjAOFAU4NCwYUP661\/\/Sl27duWNHv\/yl79EJ9fGFnnHHXdwPrwGcgIw1kTdRCMw+Crp\/vvv54m9kpwpIF46082PR4mXfnTF+TWJn861C8uRAjBpdhrbseM10fTp010\/09NPP01nnXUWPfLII66XLQXmVEC8NKdWiJfmeIk7ET\/N8lP3bgRgdJVymA\/rvmBS7ujRo6PrtzgsKu4wLNyE11aLFy\/mSb2S0q+AeJl+jb06g3jpldLenEf89EZnv51FAMZvjsj1iAKigCggCogCokBSBQRgkkokGUQBUUAUEAVEAVHAbwoIwPjNEbkeUUAUEAVEAVFAFEiqgABMUokkgyggCogCooAoIAr4TYHQAszhw4d5U0VsophK+uqrr3iTxQsuuID3O0pU3m+\/\/cafWm\/fvp1X361Tpw6dccYZfFps2Pjpp5\/yz1iMrly5cryzKtZ2URP2WcIeSB999BFvRWClYsWKRbcSQFnffvstL46HbQawOWSYErZvwD\/4YGns5P4T+ZKonLz8h+\/YGgJbP8APK6Huffzxx7zw4fXXXx+3ACJWXMbGoKgL2Cw0zMlPXsKHDz74gLC2U7169Xj\/MqTcvMSkUmwNgqXw4X8YYzHMdVfu3RsFQgswt912G2+QiA7IacKne7169Yoebu1LpK7Ii2X\/X3311Wg+7G2E1Xi3bNlCrVu35kYQ6ZxzzuFPrrHzdNu2bXNcFvZawqq7V1xxRdzfADXY4BFbFIwfPz76N+xijQ40TAmflr\/wwgu0bNkyhkUnKTdfsJ9VbMrLf3j60ksv8U7hvXv35sUOkax1faxy4Dn2ysKu5o899hgtXLgwegqsxqwujujkfoJ6jF+8hH6AylatWrGU1t5jeXkJ7wAwVoqtA0H1Q65bFPCbAqEEGEALdn7GSMfkyZN55MRuwuhN\/fr1uQzsRo0tBdD5oPPEIndWwq7SKB9w89RTT3GHhL2RPvvsM5owYQIfN2rUKN5i4Mknn+TVfbG5I\/JYCRtA4nPpkSNHUu3atenGG2\/kvNjcEenss8+m0qVLU6VKlfg8f\/vb33jRPFwPVgrGE2BswtPh2rVreeQIIwNYzTfVkSi7+qUjP7TGvlDoOKA1PjN3cl\/W0uWqL\/i9lfLyHyNld955Z8LOq3\/\/\/gyVgB+AaufOnRlcsOAh6knTpk352rErORYoRF0F5MQmjAKsWrWKdyz\/4x\/\/yKs8Y58sk5JfvEQsI16wEKUVkxbA5OYlQAcjazgW9RH+YmQ00YrZYfDSpHop9+IvBUIJMNiAERuEIeHVT+xKtgCJ1atXJ3QJIxrY0wgJ4IKRD3Q6M2bMiO4aDTDq2bNn9HiMBKCTsp7Axo4dS8899xzvXj1ixAhuFP\/973\/Tr7\/+yh0RYAPQYSVrQ0isJYOneTSE2OQRT3gY6WnQoAEfh1dKeLWEV1AYecDf8ZrC2gk79oYwuoN7RseIkQJ0mpMmTfJXzXRwNbNnz6YnnngieuS7777LUGcldPrQWU2XX345XXzxxdFfQ49kvuTlP\/zG0zkAEdAZ+\/Q9ZswYhmbUgZ9++omvd+DAgQwwqE8YZXv88cf5Hzb+xJN\/LJzgvHiFAd8s\/wA8HTp0cKCYfw\/xi5eIZcQq1nHCw8XEiROjIzC5eQmAxXHwHrGG2ELsYpQUHlspLF76t5bJlQVdgVACDEzL7RXSf\/7zHx7hSJRiO0Qs7Y8nLLwCwmiHNZzcrl07BhMr4WmtX79+BHBBJ4XGEKMnGH3BEz1W0rUACuCBkRTMl7HS3XffzU\/b+B3+pr62QD68okKnifTLL79ER1wAQ8gfOw\/GAiKM3qABxmgFnvLxSsuEBfHyeu2Anbt\/\/PHHHNaqAAAgTOaLjv8rV64kvC6MBRjMiYjdeBMQAmAGGKPe4BWhlRK9drA6dgBqx44dacqUKXytAB7Tkh+8xCrXeEgAZAIahw0bFgWY3LxEDOLBBvWtW7duPOq6Z8+eaBtg+RQmL02rm3I\/\/lBAAEaZA4MRGLxeSZTw9GS9krBA5\/bbb6esrCx+2gagADgAKFbCaAqextA5AW6mTZtGw4cP5\/kqyIcRAQAEXklgV2vMs8CoCZK1EzXgBJCCBFD64osveCIwOj80rnjiw+gLQAUjMYAWvMrCpGFruNu6HgyH4wkeDSoSzodRo5YtW\/qjRqZ4FXl1evAV\/qrp\/\/7v\/3jkykroeJAdyV8AACAASURBVPLyBfl0\/E8EMNg7CyCM1wp4fQBgAeTgtQNeU1x55ZX8M14Dwj\/ATZkyZaLXhtG6W2+9NfozRmweffRRfg1oWvKDl4inpUuX8nwm7DiPuMKrRXgF3RN5iYcB69pjPVFfL4fJS9PqptyPPxQQgFEABsO6AIxECSBidSaHDh3i+Sbo7AAJ1qsiNGqYw2IlqxPDayv8s56yMfcFQ8z\/\/Oc\/GUr+97\/\/8fAyRkbmzJnDh1vv2GPnseA9+t69exlcMjIyeCTpX\/\/6F58fE3gbNWpEgCq8IsP51Fda6MDR+WIODsqyXqVhjk3lypX9UStTuIq8Oj2MtGCDODUB3gCEVsLrurx8QT4d\/xMBDEZ3kACt+MoGmgNc77rrLgZaPOHjFSGgGP\/wuqlx48bRa8P14\/UWXi198skn\/OoSr0EBOqgPJiU\/eInRUYyAqgnx0qZNm4ReAmoQZ4hptCcY5cQrpPfee4\/+8Ic\/hNJLk+ql3It\/FAg1wKDjxysdPMVaCfMOMDk3UVIbIIyo4KkeT9VLlizhJ2aMfuDpHXMb0CmhcwToYLQEw\/5oyPD\/n3\/+OT\/NATow5+LYsWMMIbGvgwApKBOfRVuvgaxXUtjhunz58jRkyBDeawmNPb5OQmeGYW9cKxre119\/Pe7p3BrVQdkoA+UBYgBN1qRg\/1RP+1eCp+Pnn3+eh\/z\/\/Oc\/U5EiRaKFAFKskafYkq3RMet31tdcqi8oD37i1Rxev+Xmf\/Xq1bmoRAADuEXnB5\/xqgFzYVAvMHoHUEHZ8AXwAu8BmiVLloxeLiaMA3IwyRf1FJCLfBgBNOEVYKwvfvCyePHiDKtIeJDASAwePDACg9d2ibzEqz28JsQDCebNoC1APUR8x0JmmLy0H8lyhCiQXIHQAozVSQEmnH5KjZETvDKy5lXcf\/\/9PN8Fr4Dw2scadVHnreDczZo149EAHGOt+YJJwoAazIewns4BP7Gf1uL33bt3j36iiXkxaFQxEdUadYHtABlM7Ey0UzWe9AFuSLh\/vL6Ife2VvNr4N4cFDbhCdRKv7lXn5gvWdUGnZI265ea\/BZsASNQP9TNqdIDWUz0gCbAFHzEvCvMiUJ\/wag+QEvs1E64fr53we0zmRsJxqGeJPrvXvV+\/5vOTl9AIk\/XxwBD7GXVuXlp5rVgErFijb5beYfLSr3VMrivYCoQWYGAbnqww9ySVJ1d0ahi1wcJW1hdKiaoEwOO\/\/\/0vz7VQh\/oxKoBF17AWiG5C44eODqMwZ555ZvQwzHHB9ai\/V8vF0z\/Oi47StE9woTVGweBrKovZ6fii67+qPyaA4qlc\/cwb\/qGewL+8Eq4N8zNiv57SrTtByhdkL\/EaCT5hQcK8FpQMi5dBqndyrcFQINQAEwyL5CpFAVFAFBAFRAFRQFVAAEbqhCggCogCooAoIAoETgEBGA3L8FWINedAI3u+Z4n9kimvi0l0X7rH5vtNOryAoHmJ29TxRLx0WCE8Psypl7r1wOPbkdOJAvmqgACMhvxYYG7H9u1UPDOTsg8c0DhCL0u6yrukYkXCQmvJEu5r27Z\/xmW79NKqWscmK9uvfw+Kl9DPqh86fobZy1itCikV7yzl59NboCauodbxRTIz6ciBA6SbP7Y0da1nlJGqlyjf9Nj0a5sh1+VfBQRgNLyxOj2NrL7IotPhRRpEAEz8on2XXlozFADjC6M0L0LHT\/EyIqZbAGNZ4xbAWOU59TISr2bHpmY4SDZRIKqAAIxGZTAbYD5WRmDqCcBo1Akvs+h3euKlOQAT72UEYMyOTS9jSs5lhgICMBo+mg0wSxWAaSQAo1EnvMyiDzDipTkAE+9lBGDMjk0vY0rOZYYCAjAaPpoNMG8pAPMnARiNOuFlFn2AES\/NAZh4LyMAY3ZsehlTci4zFBCA0fDRbICJ7LtkpUsvbS8Ao1EnvMyiDzDipTkAE+9lBGDMjk0vY0rOZYYCAjAaPpoNMC8qANNRAEajTniZRR9gxEtzACbeywjAmB2bXsaUnMsMBQRgNHw0G2AmKwDTVQBGo054mUUfYPS8xMah2FQUO1pjb6dOnTpRZmYmrVmzhqZPn87bG2C\/LSyB\/\/3339Po0aNp\/\/79vHkkNirMa1l8L3VJFJfmAEy8lxGAMTs2vaw7ci4zFBCA0fDRbIAZpwDMYwIwGnXCyyz6AKPnJTaMXLFiBW8yiU09sSEkNoPExpJDhw7l3a2XLl3Km4h269aN98vCbtmPPfYY736NTUf9kMwGmHgvIwBjdmz6oU7JNQRLAQEYDb\/MBpihCsAMEoDRqBNeZtEHGD0vscszNhUFmIwaNYpHVDASM3fuXJo2bRpvcIpd0NeuXUtvvvkmNW\/enAoVKkStW7fmkZmbb77Zy9vP9VxmA0y8lxGAMTs2fVGp5CICpYAADBH99NNP\/MS5e\/duatmyJTfssclsgOmv3OuIQAOMaV7CHH2AsedlnTp1eLfkGTNmsOeo\/xhhQcJozPjx4+nyyy+nzZs3U5s2bfj37733HpUtW9azRi4vP80GmHgvIwAT7Nj0rNLIiUKjQOgB5sSJE\/xkWbx4capWrRo9++yz9Nlnn1HJkiWjlcBsgHlMAZhxgQUYE720BzCnvczK+oTwL9GWErt27eL6nZGRQYsWLaIJEyZQr169aPXq1TRs2DCuDw0aNKBXX32Vfv75Zwb648eP0\/Dhw3m0ZvDgwZ40kMn8NBtg4uMyAjDBjU1PKoycJHQKhB5gMPKCSYrvv\/8+m\/\/OO+\/QFVdcwcvsW8lsgOmqAMzkwAKMiV7aAxg9L\/v27csbRLZo0YLWrVtHY8eO5RjABN0lS5bQzp07qUOHDgw0mLg7Z84cuuiii2jSpEn8eunhhx\/2pKFM5qfZABPvZQRgghubnlQYOUnoFAg9wEyZMoU++OADfsePp0yMxmBCY2wyG2DuUe418johiMlEL+0BjJ6XGzduJEAM0smTJ2ngwIFUu3Zt6tevHwPN3r17GWoaNWrEE3nxxVLRokX56yS8VipTpown1SOZn2YDTLyXEYAJbmx6UmHkJKFTIPQAM2DAAJ68iK8x8Hlpnz59CF9pVKpUKW4E5mB2tqs7UaerpmHX22LFi2tBSGQDQBXW5mkdm67rT6Vc07yEFrp+OvFy3759ca9KcT78DqBSpEiRqBV4lZOdnU0lSpRIxR7bxybz0xoljY1NP39GnYqXEYAJbmzaNl8OEAU0FAg9wOArjB9++IGfMpHwZHrJJZfQgw8+GAcwO7Zv15DTH1l0Jn1GGkQATHNlBObNwAKMiV7aG4Exx0vcdzI\/zR6BifcyEq\/BjU1\/tIxyFaYpEHqAwTv\/kSNH8mTGM888k+cFPPPMM3T11VeHBGBuVQDmvcACjIle2gMYc7zEfSfz02yAifcyAjA5Y9POooRWoGMBQ3wujzlNkkSBICsQeoD57bffaMyYMTR\/\/nyeB4M1Lx5\/\/PG41UbNngNzvQIwqwMLMCZ6aQ9gzPES953MT7MBJt7LCMDkjE07ixKijF9++YWaNWvGrwnxwYIkUSDICoQeYCzzDh06xCMwse\/+rb+ZDTA1FID5R2ABxkQv7QGMeV7i\/nOLTbMBJt7LCMDkjE07ixKibcNqylguYt68eQIwQe655dpZAQEYjYpgNsD8PwVgNgUeYPKyNGhe2gMY8dLPk3hT8TICMLnHps6ihJs2baIvvviC5\/dh\/ysZgdFo\/CWLrxUQgNGwJ2idnr1JvJcpAPOtAIxGnfAyi46fkQnZ4qU5ABPvZVbWAcrK+ilHbOouSoi1fPCpPLaCwFdlGIHBqsvY40qSKBBUBQRgNJwzG2DKKQCzWwBGo054mUUfYMRLcwAm3svICEzO2NRdlPDjjz\/mbSCQ8En8c889xx8r3HDDDV5WZTmXKOCqAgIwGnKaDTDFFIA5mANgnHzpgMXR7rvvvhz7SmFfHWzXYCUMZcd+8aVhR0pZgualvdcOyb1MSTyfHWz2HJh4LyMAkzM27SxKaNmH\/a+w6rK8QvJZhZbLsa2AAIyGZEHr9HSe2CMNIl47FFYA5mgOgLHzpQMa1KysLFq+fDktXryYKleuHFc+vvb69ttv+UsIpPLly\/Mqr16loHlpD2CSe+mVzl6cx2yAifcyEq85Y9PSWXdRQi98kXOIAl4pIACjoXTQOj17AHNKAZgCOQDGzpcOhw8fph07dvCn6FhnQgUYLFGPTQQvvvhiHr72enXXoHlpD2CSe6lR3QOTxWyAifcyAjA5YzMwZsmFigJpUEAARkPUoHV69gDmiAIwRXKdA6PzpcPll1\/O5WFBQCwQqAJMt27d6MiRI3TVVVfxlg1Tp07lzzq9SkHz0h7A6Hvpld7pPI\/ZABPvZQRgco\/NdOosZYsCflVAAEbDmaB1evYAZl9Ugayscygr62zHXzosWLCAihWLvLvPDWAwnwYLBhYoUID3n9qyZQsNHz5cwwV3sgTNS3sAc9rLSIdXMnQTss2ZxBvvZRj8dCfCpZQwKSAAo+F20Do9ewCzUxmBuThHp6f7pcOaNWuiZcUCzNGjR+ngwYNUqlQpwgjMAw88QFWqVOEVkM8\/\/3yeUOhVCpqX9gAmuZde6ezFecwegYn3MgIwOWPTC53lHKKAXxUQgNFwJmidni2A+Xf8JpWX\/qFiDoBx8qVDLMB8+eWX1LNnT57Yu3TpUt6kDxN3y5UrR4MHD86xI7KGJY6z6Hh5jlL64SRnS3d+HT95QraGl46F8+GBuOcflE1WVS9KK9f9o\/LzmcrP5ys\/\/6T8\/Jvy8wUJdNmTR\/25sGLO+FKLSOQlA0yC2PShLXJJooBnCgjAaEit0+lpFONZFp0OL\/JEdylt+5cCMFfm3sC69aUDFtLCZF\/rdZNnwvx+z8l2Fk83kNgtX8dPu156qXm6zmU0wChxyfGaR2ymS2MpVxTwswICMBruGA0wXygAc03yJ0QNyXybRcdLu4CR7vzaABNCL40dgVG8ZIAxPDZ922jIhflWAQEYDWt0Oj2NYjzLotPhRUdgPlMA5joBmHQDid3ydfzkEZgQemkswChecrwaHpueNZByImMUEIDRsNJogFmpAExdARi7gJHu\/NoAE0IvjQUYxUsGGMNjU6MpliyiQJwCAjAaFcJogPlAAZibBWDSDSR2y9cGmBB6aSzAKF4ywBgemxpNsWQRBQRg7NYBowHmHQVgmgrA2AWMdOfXBpgQemkswCheMsAYHpt222XJLwrICIxGHTAaYN5UAKa5AEy6gcRu+doAE0IvjQUYxUsGGMNjU6MpliyigIzA2K0DRgPMfAVgWgvA2AWMdOfXBpgQemkswCheMsAYHpt222XJLwrICIxGHTAaYOYoANNeACbdQGK3fG2ACaGXxgKM4iUDjOGxqdEUSxZRQEZg7NYBowFmhgIw9wjA2AWMdOfXBpgQemkswCheMsAYHpt222XJLwoYPwKDPXi++OIL2rt3L6\/8WqlSJSpfvjwVLFgwh\/tLliyhwoUL08033xz3N6MB5nkFYLr4F2C88jLdQGK3fG2ACZCXCLBU\/TR6JV7FSwYYH8emdKWiQH4oYCzAfPvttzR27Fh69913WdfLLruMTp48Sdt\/3zulTZs21L17dypbtiz\/HbsiN2nShH\/3yCOPhAdgJikA081\/AOO1l3YBI935tQEmAF4isNzy02iAUbxkgEkQm9jdfdy4cbR+\/XqqUaMGderUiTIzMwkbq06fPp3OPfdcbtMqVKhA33\/\/PY0ePZr2799P9evX501UsSu8JFEgqAoYCTDz58+niRMnUpcuXahu3boMKdaIC\/bh2bNnD2G05eWXX6aPP\/6Yjhw5Qq1ateINBmvVqhUugBmnAMxj\/gKY\/PAy3UBit3xtgPG5l2gk3fTTaIBRvGSASRCbs2fPphUrVlDv3r1p5syZVLp0aWrbti01bdqUhg4dyg9s2EB14cKFvBN8xYoV6Y477qDHHnuM+vfvTzVr1gxq3yXXLQqQkQCDoL344osTviaK9fzo0aP8yujJJ5+kEiVKEOAGKdEIzMHsbMo+cMD3VaZ4ZiYVK148x47SiS6cl58fpQBMH38BTH54aRcw0plf1087Xtp5av\/ll1\/o6aefpm+++YYaNWpEd955J4O+0+Smn7hnpMPZ2XTo99j0827U52Zm0jkasZnISwaYBLG5cuVK3tUdYIJd3jGigpGYuXPn0rRp0+jUqVNUtWpVWrt2Lb355pvUvHlzKlSoELVu3ZpHZtTX5U59leNEgfxQwEiAiRUSQ6l4Qvn555\/j9B08eDAVKVKEli1bRs8\/\/zwH\/KRJk3IFmGQ7GOeHebmdU+eJnRtEAMwQBWAG+wtg8sPLdAIJ7sdu+Tp+2vHSzlP7U089xXNVMJqZlZXFnWHnzp1dqe6pxqbRIzBKXHK85hGbderU4ZHlGTNm8MPL7t27eYQFCaMx48ePp8svv5w2b95MeH2O9N5770VfobtiqBQiCnisgNEAg6e9Bg0a0O23385PKbGpR48ePPqC4Va8E77wwgvpu+++4yz33Xcf\/7OS0ZN4ByoAM8yfAOOll3YBI935tQFG00s7T+0vvPACP60jfvCKAnPFhg8fnnIz5YafRgOM4mXWx5mU9XHOkdVdu3ZRyZIlKSMjgxYtWkQTJkygXr160erVq2nYsGHsE9rAV199lR\/iMFJz\/Phx9hCjNXiQkyQKBFUBowFm1apVNGTIEHr\/\/fdz9QcNAF4lIaGBxpBr165d+V1yKADmcQVgnvEnwHjpZbqBxG752gAT42XW6kzKWp33q0Tdp3bEAeZRDBgwgF566SWqUqVKyu2dG34aDTBKXPIITILY7Nu3L8\/ba9GiBa1bt44\/XMBEXUzQxTy\/nTt3UocOHRhoMHF3zpw5dNFFF\/FoM9q6hx9+OGUvpQBRIL8UMBpgMPTdsGFDeuKJJzh4zzrrrKjOGH1RE4Ia82BC9RXSowrAjPcnwHjppV3ASHd+bYDR9FL3qX3BggX8ZI8JoocOHeKndTzBu5Hc8NNogFG8ZIBJEJsbN24kQAwSvrIcOHAg1a5dm\/r168dAg+UjADWYv4SJvPhiCXOY8HUSXiuVKVPGDTulDFEgXxQwGmAAI5h0uGHDhhzi4nfnnXeeluhGv0J6UAGYKf4EGC+9TDeQ2C1fG2A0vdR9asccFeTFaKQK9VqBk0cmN\/w0GmAULxlg8ojNffv28auk2ITfAVQw189K0D07O5s\/WpAkCgRdAaMBBmsj4PPoefPmxb0SgmkYRtVdA8FogLlfAZip\/gQYL720Cxjpzq8NMJpe2nlqv\/baa+nXX3+Njl5iPgW+Sko1ueGn0QCjeMkA49PYTLUuyPGigFMFjAYYTBTEjPtPPvkk6SfVeQmYmTmFDmZfrKVxAbo9R75T9LbWscjkxvE6HR43iPgK6R4FYGb4E2Dc8hL3\/N32iXn6cebvfy3yu5eHY3In8tLKb2X77ff\/SeQl\/nSGUh+s\/OpFWcfr+OnES92ndu3KayOjG37inv+jeKl6EfuiuALdTj8q13hA8aKQ8vfjMT+XoNtJ9eqXBLEdmVF3OlnHFKbb6cKKyeMrkZccrz6NTRu2S1ZRwFUFjAaYAwcO8LoHeOfbuHHjuCFWLOaE9\/s66aKLmtD3u3frZCWirxLku0rzWGRL\/XidDi8KMH9WAGZu8gbWxs24ltUtLyOjaefmeV3WgPsZv3sRCzBEOb08PUAfKfZItPREXhIVUco4nV+9rMjxOn5ypxcQL3FPbviJe\/5e8VL1IvZFyVn0VQ6AOal4UUyx4GDMz4XoqxwAc06C+rBfKcPyF\/VJG2AULzlefRqbrgW5FCQK2FTAaIDBDHxM4E2UpkyZEvduOC\/djH6F1FoBmPk5AcbOwmeWjphMiE\/R3Zr06aWX6X4lZLd8bYDR8NJm+5C27G74afQrJMVLBpgEsZk2g6RgUSAAChgNMG7pbzTA3KYAzKKcAGNn4TPMr8CCZ8uXL6fFixdT5cqV3bLBlXJ0vLQLGOnOrw0wGl66IqJPCjEaYBQvGWASxKZPrJDLEAXyRQGjAQafFT733HO8HgJm3scmrGuBGfo6SafT0ynHqzw6HR43iHjtcKsCMO\/lBBg7C58dPnyYduzYQY8\/\/jivNeEWwHjpZbqBxG75On7qeulVHUx2Hjf8NBpglLjkeE0Qm8l0lr+LAiYrYDTAfPXVVzwHZsSIETnWO7j++uu1J\/YaDTA3KQDzYe5zYOwsfIaFtUaOHOkawHjppV3ASHd+bYCx4WV+N2pu+Gk0wCheMsDkEZv57aecXxTIDwWMBhgs5NSzZ0\/ecTqVZDTA3HAaYLL+k0lZ\/3G+XDkWPitWLDIN0m2A8dLLdAOJ3fK1ASbGS+7wVvlzQjauzQ0\/jQYYxUu\/+5lK+yrHigJOFTAaYLBFQLt27ejKK6\/k\/UBiV+LFNvIFCxbU0i1zyhQ6eLHmZ9S3J\/iM+m0bn1G7cLxOhxd9hfRHZQTm85ydnp2FzyxB3QYYt7zkz6gnan5G\/bsXcZ9RJ\/Ay18+oE3gJfc5Qysj1M+rfj9fxk18haXipVeE9yOSGn\/wZteJlnp9R357gM2rFizw\/o749wWfUCepDrp9R327jM2rFS47XBLHpgVVyClHAtwoYDTCYj3HjjTfSjz\/+mGMhOztzYC5q0oS+P\/\/8XE2M\/cs5c+Zwvp2xudu3z3GsWtpPVo7fj487oH17yjW\/WvKcOVqf3UYBppoCMBtyAoydhc\/SBTBuecmjabVrx6l2tqLhNb\/\/XPF3L2Lx81gCL638VjFf\/P4\/GYm8JKIqShlWfuv4\/yl1QRtgNLz0S0vkhp\/w8hfFy7LKDd4a83P9OXPoZeXv\/1G8qKf8PXbstnwCP8snqA\/vKWX89\/ef0TacobsOjOIlx2uC2PSLn3IdokB+KGA0wHz66ae8WRk2jku095Gu4MleISUam4kDmAQnUo9xM79OhxcFmMoKwGzO\/bVDfi58lk4viyv+qGNoycbP3M4fP93cxjowNrzUrfvpyueGn4jLk9vj6+8VygWrWxWqAKPe373KL6YnEeCuBH9\/Vvnd1zE\/awOM4iXHax6xmS6fpFxRwM8KGA0wP\/30E2\/i+MEHH+TYJ8SOKUYDzP8pAPONP+dNpNNLYwAmIF4i9tzw02iAUbxkgPFpbNppSyWvKOCmAkYDzO7du3kzxz179hDmvJQte3qAefjw4a4tZBfoEZiLFIDZ5U+ASaeXxgBMQLxEA+aGn0YDjOIlA4xPY9PNDknKEgXsKGA0wOzfv5\/wZUyidM8992i\/VjJ6BOZCBWB+8CfApNNLYwAmIF4iHt3w02iAUbxkgPFpbNrpcCSvKOCmAkYCzIcffkiXXXYZlS9fPk+tsCPu1VdfnVRPowEmUwGYA\/4CGC+8NAZgfO4lAs1NP40GGMVLBhifxWbShlMyiAJpVsBIgPnss8+od+\/e9Mc\/\/pFuuukmuuiii\/gf9vTBKrHYCfett97in994442kEhsNMOcqAHPIXwDjhZfGAIzPvUSguemn0QCjeMkA47PYTNpwSgZRIM0KGAkw0OzIkSP02muv0bJlywgjLfhs00p169al1q1bU7NmzeiMM85IKrHRAFNQAZgT\/gIYL7w0BmAC4KWbfhoNMIqXDDAJYtPORqu\/\/PILPf300\/TNN99Qo0aNeH5g0aJFk7Z\/kkEU8KsCxgKMKvjevXsZYjCR1+4n1UYDzCkFYAr4D2DS7aUxABNAL+Gt09g0GmAULxlgEsSmnY1Wn3rqKTp48CB16dKFN1ytWrUqde7c2a99k1yXKJBUgdAATFIl8shgNMAcUwAmw\/8A47aXxgBMCL00dh0YxUsGmASxaWej1RdeeIFHnsuVK0czZ86kLVu2EL7GlCQKBFUBARgi2rp1Ky1fvpzOO+88aty4MZUoUSLOT6MB5pACMOcGG2CceGkMwBjmJYIwLz+NHoFRvMzKyKSsjJz7lFkNlZ2NVrEK+YABA+ill16iKlWqBLXvkusWBSj0ALNz505q2rQpL3iH10tTp06l1atXx+1ebTTAHFAAJjO4AOPUS2MAxiAv0TYn89NogFG85BGYBLG5a9cuXqQzIyODFi1aRBMmTKBevXpxGzZs2DDu4rAPHJaTQB583HDo0CEaPHgwVaxYUbpAUSDQChgNMJjIi0awUqVKcSZ99NFHdP311\/NmjvPmzaM1a9bwO2Gku+++mxo2bMj\/tZLRAPO9AjBl\/Akw6fTSGIAJiJeIKzf8NBpgFC8ZYBLEpp2NVpG3dOnS9MgjjwS605KLFwUsBYwGmM2bN\/OXRnfddRc98cQT0d2n0fBt2LCBXxlhQa1Tp07xUwz+H59d4\/1w7PowyH8wO5uyDxxIWHP8tBJv8cxMKlY896Hm2BvgHYz\/owBMeX8CTDq99DPA6PoZJC9RB93wE\/eMdCo7m079Hpt+3gupQGYmFdCIzUReMsAkiE07G61ee+219Ouvv9JZZ50VHZnBV0mSRIGgKmA8wLRp04YyMzOpQoUK9Oyzz1KxYsUoFmAs495++20aOnQotWrVivCkonb0O5RN42L\/7ieAwXXZ2czx38p9\/UFjt9z8qOzo8NLlpZ8BRtdP1OmgeGkBTKp+mjwCo3oJzfKKzfzcaDU\/2gM5pygABYwHmG7dutH8+fPpoYceIrwvxsS1W265JToCc+LECX4fjBEZvDOuXr16jpph8iukfykAc6WPASZdXpoCMEHx0gKYVP00GWBUL6GZX2NTulJRIL8UCAXA4Auj48eP8yeDr7\/+Oq8HY71CwuS2F198kX9fqFAh9qFAgQL8z0omA8wXCsBc43OASYeXpgBMULyMBZhU\/DQZYFQvoZlfYzO\/Oi85ryhgNMBgx9sZM2ZQ\/\/79o05j4Se898UsfaxCib+98sorcTUBs\/i7du0aCoD5TAGY63wKMOn00hSACYqXCCw3\/DQZYFQvoZlfY1O6UVEgvxQwGmDcEtXkZwKtmgAAE59JREFUEZiVCsDU9SnApNNLUwAmjF6aupCd6iXqv+mx6VaMSznhUUAARsNrkwHmAwVgbhaAoduVOvF2kjridv5s5Xw6k7JRR3W9tLN\/jnUpAwcOpPvuu89Xa4eYPAKjegkfTI9NjaZYsogCcQoIwGhUCJMB5h0FYJoKwAQWYHS9tLN\/Dj7TxRpJmKuyePFiqly5skbEeJPFZIBRvYSipsemN7VGzmKSAgIwGm6aDDBvKgDTXAAmsACj66Wd\/XMw4X3Hjh30+OOP06RJkwRglPbirgTtx7PK776O+fkMjfhCe6N6iSJMj02NpliyiAIyAmO3DpgMMPMVgGmt0cDa1c9P+RN5acocmFgvX83MpFeTLJpmZ\/+cFi1a0MiRIwVgPAIYNS5xWtNj00\/thFxLMBSQERgNn0wGmDkKwLQXgAnsCIyul3b2z8HCj0gCMIkbinSNwKhe4uymx6ZGUyxZRAEZgbFbB0wGmBcVgOmYAGDsTPrEHjcTJ06kTZs25dhTCrpjRV2siGylTp06xW3bYNcbu\/lNHoHR8RJ62dk\/x9JXAMZbgFG9xNkTxabd+i\/5RQGTFJARGA03TQaY5xWA6ZIAYOxM+hw7dixt3bqVOnbsSPhyZdCgQVS3bt2oylgV+dtvv+U9qpDKly\/P6\/F4lUwGGB0vobOd\/XMEYPKumekagVG9xFUkik2v4kbOIwr4UQEBGA1XTAaYSQrAdEsAMHYmfWIvqREjRvCWDFhEMDs7m3r06BFVGYCTkZFBF198Md1www1UokQJDQfcy2IywOh4Gatk0PfPMfkrJNVL+JYoNt2LDClJFAieAgIwGp6ZDDDjYgDmvcxMei+PiZ86kz6bNGlC69ev500zly1bRgsXLuSvV6yE\/W\/wmumqq64ijOxMnTqVqlWrpuGCO1lMBphYL6HWYyGYz2TqQnaql2Hw050Il1LCpIAAjIbbaqen7j793alTOUqpELOXUqJTqMekmn9nzEl0Fj5DdtzXKGUEpk+CTs\/OpM+WLVvSrFmzqEyZMrRkyRLecyp2d2\/Mp8GeU9hraubMmbRlyxbeo8qrhHs+qNyzuhDddMXPB5J4+bcU89+rlK8unFdMA0R0vfRKZy\/Og3uupHj5sHLiZqeOxP3m6wJF8ry0KxQvv07i\/RVK+Sh8sXKO2M+qtzr0EuUmik0vdJZziAJ+VUAARsMZkwFmiNIBDE7QwNqZ9Im89erVo8aNG\/M+UzVr1uT5LgcPHqRSpUoRRmAeeOABqlKlCo0ZM4bOP\/98ni\/jVTIZYHS89EpnL85jMsCoXkLPRLHphc5yDlHArwoIwGg4YzLADFQAZlgCgLEz6RNfGd199908MbdkyZI8Dwa\/69mzJ6\/munTpUho1ahT\/vVy5cjR48GDO51UyGWB0vPRKZy\/OYzLAqF5Cz0SxaecLQcsTP24L4UV9kXOYp4AAjIanJgPM4wrAPJPHELfupM\/jx4\/T\/v37+TVSonTixAnCCq\/WGiMaFriWxWSAseOla4LmY0EmA4zqJWROFJt2vhD087YQ+ViN5NQBVkAARsM8kwHmUQVgxmu8o9eQzLdZTAaYMHpp6hwY1UsEVKLYtPOFoJ+3hfBtgyEX5msFBGA07DEZYB5UAGaKAAwFdRJvGL00FWBUL9FM5RWbOl8IXn755dza+XFRQo1mWLKIAjkUEIDRqBQmA8z9CsBMFYAJLMCE0UtTAUb18ovMTPoiwRIHdr4Q9PO2EBrNsGQRBQRgnNQBkwHmHgVgZgjABBZgwuilqQCjeol2K1Fs2vlC0Gr7ZATGSS8gx\/hRARmB0XDFZID5swIwcwVgAgswYfTSVIBRvUQzlSg27XwhKACj0dhLlkApIACjYZfJAHOHAjBvCMAEFmDC6KWpAKN6iWYqr9jU\/UJQo7mTLKJAYBQQgCGi3bt387L3BQsWpFtvvTXH\/jwmA8xtCsAsCjjA6Hhp6kq8pnmJVjQvP03+jFr1EloEPTYD0yvKhQZGgdADDFaIbdiwIV1zzTUMLmvXrqXFixfzhoNWMhlgblUA5r0AA4yul6YCjEleIvaS+WkywKheQo8gx2ZgekS50EApEHqAwX48GH3Bf5GwDP7kyZN5qfswAMxNCsB8GGCA0fXSVIAxyUvEXjI\/TQYY1UvoEeTYDFSvKBcbGAVCDzBYyv7ss8+mPn36sGmdO3emRo0aUevWrUMBMLUVgPkkwACj66WpAGOSlwi+ZH6aDDCql9AjyLEZmB5RLjRQCoQeYLC5IDYcxP49SAMGDKBLLrmEQSZ2BOZgdjZlHzjAv\/LzbtTFMzOpWIL1IhLVSnQAf1QA5vMAA4yul0ezs+nX372ELn7ejfqszEwqrOGnaV7Cl2R+4p6RLsjOpgt+99PPu1H\/mJlJPzr0EvcZ5NgMVK8oFxsYBUIPMFOmTKGff\/45OgLTpUsXateuHd10002hGIGppgDMhgADjK6Xpo7AmOQlgi+ZnyaPwKheQo8gx2ZgekS50EApEHqAWbVqFQ0ZMoQn7mJVy+bNm\/OOybEbEZo8ibeyAjCbAwwwul6aCjAmeYlWNJmfJgOM6iX0CHJsBqpXlIsNjAKhBxjsnIyh6jVr1vAOyf369aNOnTrFGWgywPyfAjDfBBhgdL00FWBM8hIBmMxPkwFG9RJ6BDk2A9MjyoUGSoHQA4zl1t69e3ky77nnnpvDQJMB5iIFYHYFGGB0vTQVYEz0Ep7mFpsmA4zqJXQwITYD1TvKxfpeAQEYDYtMBpgLFYD5wQCAyctSeGkqwITRS1NX4lW9RJ02PTY1mmLJIgrEKSAAo1EhTAaYTAVgDgjABHYrgTB6aSrAqF6imTI9NjWaYskiCgjA2K0DJgPMuQrAHBKACSzAhNFLUwFG9RJtlumxabddlvyigIzAaNQBkwGmoAIwJwRgAgswYfTSVIBRvUQzZXpsajTFkkUUkBEYu3XAZIA5pQBMAQGYwAJMGL00FWBUL9FmmR6bdttlyS8KyAiMRh0wGWCOKQCTIQATWIAJo5emAozqJZop02NToymWLKKAjMDYrQMmA8whBWDOTQAwx44do3HjxtH69eupRo0avE5OZmYmr50zffp0\/vS8e\/fuVKFCBTpy5AhNnDiRNm3axLt8W1s02NU8XflN\/gpJx0voaoqfJn9GrXoJ3xLFZrriRMoVBYKggIzAaLhkMsDsVwCmRAKAmT17Nq1YsYJ69+7NOwSXLl2a2rZtS02bNqWhQ4fS9u3befXihQsX0tixY2nr1q3UsWNHGjhwIA0aNIjq1q2robI3WUwGGB0vobIpfpoMMKqX8C1RbHoTNXIWUcCfCgjAaPhiAQw2SsSGjm5t5piVlUU9evTgK6hQoECeV\/LdqVNxf1fz7yQi6\/ou0XwNhPv6XgGYMgmOXblyJZUrV44qVqxIo0aNogIFCvBIzNy5c2natGl06tQpqlq1Kq1du5ZatWpFI0aMoOrVq9OMGTMoOzs7eo8aUqc9iwUw2CTR2tDRjc0cY718IImXf1O8vFfJ\/\/bvKljXWEzDT10vUbQpfsYCDDZKxIaObm3maPn5dRIvrzh1JEedXVygSNzvniUi6\/q2OvQSBSaKzbQHjJxAFPCxAgIwGua0b9+ePv30U42c\/sjy66+\/0n\/\/+9+kF5PovmrVqkVz5sxJeGydOnVoz549DCbbtm2j3bt3U\/\/+\/TkvRmPGjx9PTZo04VdNxYoVo2XLlvGozKRJk5Jei1cZguYldMnLE0u3RPd19NdfaXce9SDofobJS9164FUcyXlEAT8oIADjBxd8fg3Y5LJkyZKUkZFBixYtogkTJlCvXr1o9erVNGzYML76Bg0a0IIFC6hly5Y0a9Ys3gxzyZIltGHDBurbt6\/P7zBclyd+hstvuVtRwFQFBGBMddbF+wKAYBSgRYsWtG7dOp7nMnr0aJ7nAkjZuXMndejQgSf1Im+9evWocePGPDpTs2ZNhhpJ\/lFA\/PSPF3IlooAo4FwBARjn2oXmyI0bN0ZHUU6ePMmTc2vXrs07dwNosNkeoKZRo0a0efNm\/vKoaNGiPGqD102FCxcOjVZBuFHxMwguyTWKAqJAMgUEYJIpJH+PKrBv3z6GktiE3+Ez6iJFTk9cPH78OO3fv59fI0nyrwLip3+9kSsTBUSB5AoIwCTXSHKIAqKAKCAKiAKigM8UEIDRMARf2+CLmoIFC9Ktt95KJUqU0DgqcRaMTLzzzjtxC7w5LR\/rrSxfvpzOO+88nnNiXZfT8hzfVMAOdEsf8TL\/jfejl1BFYjP\/64ZcgfkKCMAk8fjgwYO8ouw111zDgIC1ThYvXsxf5NhJWCvlyy+\/5HVT8AkyIAbJafmYOItPl+vXr09ly5alqVOn8ldBZ599tivXa+fegpTXqd6x9yhe+sNxP3oJZSQ2\/VE\/5CrMV0AAJonHWHkWoy\/4LxK+sJk8eTJVqVLFVu3A8u3dunUjzDs4evRoFGCclj9v3jz+6gcLbiFh4ixAC8mN67V1cwHK7FTv2FsUL\/1huB+9hDISm\/6oH3IV5isgAJPE48GDB\/OoRp8+fThn586d+Wub1q1bO6odgA6snWKNwDgtH68vMBKASbX4\/5tuuokhCwvHuXm9jm7Sxwc51TvRLYmX+Wu0H72EIhKb+Vsv5OzhUUAAJonXGDXBWibWpoQDBgygSy65hEHGSVI7vVTLf\/vtt3k\/Iizhj\/U9Ui3PyT0F6Rg39REv89d5P3sJZSQ287d+yNnNV0AAJonHU6ZMoZ9\/\/jk6AtOlSxdq164dj3g4SWqn57T8EydOEJ5AsdItRnSw9xCS0\/Kc3EsQj3FTH\/Eyf2uAH72EIhKb+Vsv5OzhUUAAJonXq1atoiFDhvDEXSzB3rx5c9552ekaJ2qn57R8LNv\/4osv0uuvv06FChXiu8Ami5jI6+b1mhYKTvXWeYXktGzx0lktc6p3Or1E2eKnMz\/lKFHArgICMEkUw6JsGKoGeBw+fJhXn+3UqZNdnaP5VYBxWj6W6X\/llVfirgP7E+Ha3Lxexzfq0wOd6q3T6TktW7x0Vlmc6p1OL1G2+OnMTzlKFLCrgACMpmJYLh+TY7HqbDqS2+W7XV467jk\/y0ynPm6X7XZ5+al7Os6dTn3SUXY6ykyHrlKmKOB3BQRg\/O6QXJ8oIAqIAqKAKCAK5FBAAEYqhSggCogCooAoIAoETgEBmMBZJhcsCogCooAoIAqIAgIwHtSBJUuW8H5FN9xwg2tn27hxI23atInatm3rWplSUHIFxMvkGgUlh3gZFKfkOkWBxAoIwKS5Zhw4cICaNWtG7777LhUtWtS1s2E7AqxFM3v2bKpYsaJr5UpBuSsgXppTO8RLc7yUOwmvAgIwDrzH55vr16+nChUqUKlSpXjhKutnLO0fm8aOHct7H+HTSmzyhrzYAuCbb76ha6+9ls4880zeIBJfN1133XW8lgtGVy688EJed+bQoUO8EjCWJ8eidYCVypUr8yn+\/ve\/0+bNm2ncuHEO7kIOgQLipTn1QLw0x0u5E1FARwEBGB2VEuTBirw\/\/PADzZ8\/nxeUA0xggbvMzMy43LfccgsNGjSI6tatS8888wy99tpr\/PciRYrQkSNHeFTm\/PPPp6+\/\/poeeugheuCBB6hFixa0e\/duKl68OO3Zs4e3LtixYwddeeWVtG7dOi6jRo0aDC8Y3cGrJJQnyZkC4qUz3fx4lHjpR1fkmkSB9CggAONQV+wqffPNN9Of\/vQnmjt3Lm+kqM5xwcgLRkuWL1\/OIycAmJdffpk+\/\/xzfvKvWrUqPfroo\/Twww\/TtGnT6P333+fF6QAwGMkBFK1cuZLuuecemjNnDtWqVYsXqqtWrRovVofRGZSB3W8xmiPJmQLipTPd\/HiUeOlHV+SaRIH0KCAAk4Ku2KytR48e9Oc\/\/5mGDx+eo6StW7dS48aNeXQlIyODAQavkZ599lnOe+mll\/JWAFdffTVv\/Pb888\/TokWLGGDatGlD7du359dJgCSUVbBgQcJmkpgQjI0bkQAu2MyxSZMmKdyJHCpemlMHxEtzvJQ7EQXyUkAAJoX6AZjA6Mtll13G4FG4cOG40iz4+PLLL\/lVEQAGr42wCaMFMG+99RZVqVIlB8BgKBxQYpWxbds2PkYFmKuuuoqmTp3K82ckOVdAvHSund+OFC\/95ohcjyiQHgUEYBzq+uGHH9J9\/7+9O8ZREIjCADwlB+MC9jZGGq3suIKJl+AEEAuOwWG8wWYmWbMbTUwIzxj8aHXeyPcs\/sDA7PdpGIZ0OBzKJo9t2\/6rlnexzrtEX6\/XlIPG0gHmdruVqzfTND2svZl5Wl85TC\/X03a9XE8vnQmBVwICzCuhJ5\/n4FDXddrtdul4PJaNHrfb7X1x7d8hm82mfC\/fBpoTYH4X6j67ApOfSmqapizsdcwT0Mt5bp84Si8\/sSt+E4E4AQEmzrZU7vs+5dtEXdctPtP5fE5VVaXT6bR4bQUfBfRyPf8KvVxPL53J9woIMMG9z+99yYtyL5fL\/f0tS0yZX8SVb1uN41gW9TriBfQy3vhdM+jlu6TNQyBOQICJs1WZAAECBAgQCBIQYIJglSVAgAABAgTiBASYOFuVCRAgQIAAgSABASYIVlkCBAgQIEAgTkCAibNVmQABAgQIEAgSEGCCYJUlQIAAAQIE4gQEmDhblQkQIECAAIEgAQEmCFZZAgQIECBAIE5AgImzVZkAAQIECBAIEvgBW1ODBm9w\/QUAAAAASUVORK5CYII=","height":337,"width":560}}
%---
