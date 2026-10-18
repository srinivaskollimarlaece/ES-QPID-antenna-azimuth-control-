function sqpid_figures(opt, dpi)
% SQPID_FIGURES  Draw and save every figure of the article from results_option<opt>.mat
%   sqpid_figures(opt)       saves PNG files at 300 dpi in the folder figures_option<opt>
%   sqpid_figures(opt, dpi)  chooses the resolution

  if nargin < 2, dpi = 300; end
  S = load(sprintf('results_option%d.mat', opt));  R = S.R;
  odir = sprintf('figures_option%d', opt);
  if ~exist(odir, 'dir'), mkdir(odir); end
  Ls = R.Ls;  main = R.main;  allc = R.allc;  prop = R.prop;
  nm = numel(main);

  % colour and line style per controller
  sty = struct();
  sty.PID      = {[0.90 0.45 0.00], '-.', 1.4};
  sty.S_PID    = {[0.00 0.45 0.74], '--', 1.4};
  sty.FS_PID   = {[0.20 0.63 0.17], '-',  1.2};
  sty.S_QPID   = {[0.70 0.05 0.15], '-',  2.0};
  sty.ES_QPID  = {[0.70 0.05 0.15], '-',  2.0};
  sty.QPID     = {[0.50 0.50 0.50], ':',  1.4};
  sty.ES_noweight = {[0.55 0.35 0.65], '--', 1.2};
  sty.ES_nofilt   = {[0.35 0.60 0.80], '-.', 1.2};
  sty.ES_noDw     = {[0.55 0.55 0.20], ':', 1.4};
  key = @(n) strrep(n, '-', '_');
  cgust = [1.00 0.95 0.80];
  ton = 8;  toff = 14;

  %% ---------------- step responses ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,1,j);  hold on;
    patch([ton toff toff ton], [-0.1 -0.1 1.45 1.45], cgust, 'EdgeColor', 'none');
    h = zeros(1, nm+1);
    h(1) = plot([0 24], [1 1], 'k:', 'LineWidth', 1.0);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i+1) = plot(tr.t, tr.y, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.45]);
    ylabel('Azimuth (deg)');
    title(sprintf('(%s) L = %.1f s', char('a'+j-1), Ls(j)));
    text(ton+0.15, 1.36, 'wind gust', 'FontSize', 8);
    if j == 1, legend(h, [{'Reference'} main], 'Location', 'southeast'); end
    if j == 3, xlabel('Time (s)'); end
  end
  savefig(fh, odir, 'fig_step_response', dpi);

  %% ---------------- control effort (left: start-up, full scale; right: regulation, zoomed) ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,2,2*j-1);  hold on;
    h = zeros(1, nm+1);
    h(1) = plot([0 4], [20 20], 'k:', 'LineWidth', 1.0);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i+1) = plot(tr.t, tr.u, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([0 4]);  ylim([-2 22]);
    ylabel('u (V)');
    title(sprintf('(%s) L = %.1f s, start-up', char('a'+j-1), Ls(j)));
    if j == 1, legend(h, [{'Actuator limit'} main], 'Location', 'northeast'); end
    if j == 3, xlabel('Time (s)'); end
    subplot(3,2,2*j);  hold on;
    patch([ton toff toff ton], [-0.6 -0.6 2.6 2.6], cgust, 'EdgeColor', 'none');
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      plot(tr.t, tr.u, '-', 'Color', s{1}, 'LineWidth', 0.8);
    end
    box on;  grid on;  xlim([5 24]);  ylim([-0.6 2.6]);
    title('regulation and wind gust (zoomed)');
    if j == 3, xlabel('Time (s)'); end
  end
  savefig(fh, odir, 'fig_control_effort', dpi);

  %% ---------------- pointing error around the gust ----------------
  fh = newfig(7.2, 7.4);
  for j = 1:numel(Ls)
    subplot(3,1,j);  hold on;
    yl = [-0.33 0.37];
    patch([ton toff toff ton], [yl(1) yl(1) yl(2) yl(2)], cgust, 'EdgeColor', 'none');
    plot([6 24], [0.02 0.02], 'k:', 'LineWidth', 0.8);  plot([6 24], [-0.02 -0.02], 'k:', 'LineWidth', 0.8);
    h = zeros(1, nm);
    for i = 1:nm
      ii = find(strcmp(allc, main{i}));  tr = R.traj{ii,j};  s = sty.(key(main{i}));
      h(i) = plot(tr.t, 1 - tr.y, s{2}, 'Color', s{1}, 'LineWidth', s{3});
    end
    box on;  grid on;  xlim([6 24]);  ylim(yl);
    ylabel('Pointing error (deg)');
    title(sprintf('(%s) L = %.1f s', char('a'+j-1), Ls(j)));
    if j == 1, legend(h, main, 'Location', 'northeast'); end
    if j == 3, xlabel('Time (s)'); end
  end
  savefig(fh, odir, 'fig_gust_error', dpi);

  %% ---------------- performance indices ----------------
  fh = newfig(7.2, 3.4);
  idx = [1 2 3];  nmI = {'ISE (deg^2 s)', 'IAE (deg s)', 'ITAE (deg s^2)'};
  for q = 1:3
    subplot(1,3,q);
    V = zeros(numel(Ls), nm);
    for i = 1:nm, V(:,i) = squeeze(R.nom(find(strcmp(allc, main{i})), :, idx(q))); end
    hb = bar(V, 0.85);
    for i = 1:nm, s = sty.(key(main{i}));  set(hb(i), 'FaceColor', s{1}); end
    set(gca, 'XTick', 1:3, 'XTickLabel', {'0.3', '0.8', '1.5'});
    xlabel('Delay L (s)');  ylabel(nmI{q});  grid on;  box on;
    if q == 1, legend(main, 'Location', 'northwest'); end
  end
  savefig(fh, odir, 'fig_indices', dpi);

  %% ---------------- parametric robustness ----------------
  fh = newfig(7.2, 6.6);
  for j = 1:numel(Ls)
    for q = 1:3
      subplot(3,3,(j-1)*3+q);  hold on;
      h = zeros(1, nm);
      for i = 1:nm
        s = sty.(key(main{i}));
        v = squeeze(R.rob_IAE(i,j,q,:));  v(squeeze(R.rob_ok(i,j,q,:)) < 0.5) = NaN;   % unsettled runs not drawn
        h(i) = plot(100*R.lev, v, [s{2} 'o'], 'Color', s{1}, 'LineWidth', s{3}, 'MarkerSize', 3);
      end
      box on;  grid on;  xlim([-32 32]);
      if j == 1, title(sprintf('%s error', R.par{q})); end
      if q == 1, ylabel(sprintf('IAE,  L = %.1f s', Ls(j))); end
      if j == 3, xlabel('Error (%)'); end
      if j == 1 && q == 1, legend(h, main, 'Location', 'north'); end
    end
  end
  savefig(fh, odir, 'fig_robustness', dpi);

  %% ---------------- Monte Carlo ----------------
  fh = newfig(7.2, 3.4);
  lab = {'IAE (deg s)', 'ITAE (deg s^2)', 'Gust peak error (deg)'};  col = [2 3 5];
  for q = 1:3
    subplot(1,3,q);  hold on;
    M = zeros(numel(Ls), nm);  Sd = M;
    for i = 1:nm
      for j = 1:numel(Ls)
        x = squeeze(R.mc(i,j,:,col(q)));  x = x(isfinite(x));
        M(j,i) = mean(x);  Sd(j,i) = std(x);
      end
    end
    hb = bar(M, 0.85);
    for i = 1:nm, s = sty.(key(main{i}));  set(hb(i), 'FaceColor', s{1}); end
    wgrp = 0.85;
    for i = 1:nm
      xo = (1:numel(Ls)) + (i - (nm+1)/2)*wgrp/nm;
      for j = 1:numel(Ls)
        plot([xo(j) xo(j)], [M(j,i)-Sd(j,i) M(j,i)+Sd(j,i)], 'k-', 'LineWidth', 1.0);
      end
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', {'0.3', '0.8', '1.5'});
    xlabel('Delay L (s)');  ylabel(lab{q});  grid on;  box on;
    if q == 1, legend(hb, main, 'Location', 'northwest'); end
  end
  savefig(fh, odir, 'fig_monte_carlo', dpi);

  %% ---------------- step-size study ----------------
  fh = newfig(7.2, 3.4);
  cs = [0.10 0.10 0.10; 0.70 0.05 0.15; 0.00 0.45 0.74; 0.20 0.63 0.17; 0.90 0.45 0.00];
  ls = {'-', '-', '--', '-.', ':'};
  for j = 1:numel(Ls)
    subplot(1,3,j);  hold on;
    h = zeros(1, numel(R.steps));  lg = cell(1, numel(R.steps));
    for a = 1:numel(R.steps)
      tr = R.stptraj{j,a};
      h(a) = plot(tr.t, tr.y, ls{a}, 'Color', cs(a,:), 'LineWidth', 1.3);
      lg{a} = sprintf('%g deg', R.steps(a));
    end
    box on;  grid on;  xlim([0 10]);  ylim([0 1.15]);
    xlabel('Time (s)');  title(sprintf('L = %.1f s', Ls(j)));
    if j == 1, ylabel('y / step size'); end
    if j == 3, legend(h, lg, 'Location', 'southeast'); end
  end
  savefig(fh, odir, 'fig_step_size', dpi);

  %% ---------------- Popov plot ----------------
  fh = newfig(7.2, 3.4);
  for j = 1:numel(Ls)
    s = R.pop{j};
    subplot(1,3,j);  hold on;
    X = real(s.H);  Y = s.w.*imag(s.H);
    plot(X, Y, '-', 'Color', [0.70 0.05 0.15], 'LineWidth', 1.6);
    sel = X > -1/s.k - 0.6;
    xl = [-1/s.k - 0.45, max(0.6, 1.08*max(X))];
    yl = [max(1.12*min(Y(sel)) - 0.05, -3), 0.5];
    q = s.q(1);                                  % regulation mode (first entry)
    if q > 0
      xx = [xl(1), xl(2)];  plot(xx, (xx + 1/s.k)/q, 'k--', 'LineWidth', 1.0);
    else
      plot([-1/s.k -1/s.k], yl, 'k--', 'LineWidth', 1.0);
    end
    plot(-1/s.k, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    box on;  grid on;  xlim(xl);  ylim(yl);
    xlabel('Re H(j\omega)');  title(sprintf('L = %.1f s', Ls(j)));
    if j == 1, ylabel('\omega Im H(j\omega)'); end
    if j == 3, legend({'Popov locus', 'Popov line', '-1/k'}, 'Location', 'southeast'); end
  end
  savefig(fh, odir, 'fig_popov', dpi);

  %% ---------------- Option 2 only ----------------
  if opt == 2
    % ablation of the design cost
    fh = newfig(7.2, 3.0);
    names = {'PID', 'S-PID', 'FS-PID', 'S-QPID', 'ES-noweight', 'ES-nofilt', 'ES-noDw', 'ES-QPID'};
    V = zeros(numel(Ls), numel(names));
    for i = 1:numel(names), V(:,i) = R.J(find(strcmp(allc, names{i})), :)'; end
    hb = bar(V, 0.9);
    for i = 1:numel(names), s = sty.(key(names{i}));  set(hb(i), 'FaceColor', s{1}); end
    set(hb(4), 'FaceColor', [0.95 0.60 0.65]);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'L = 0.3 s', 'L = 0.8 s', 'L = 1.5 s'});
    ylabel('Robust design cost J');  grid on;  box on;
    legend(names, 'Location', 'northwest');
    savefig(fh, odir, 'fig_ablation', dpi);

    % probability weight and its effect (L = 0.8 s)
    fh = newfig(7.2, 4.6);
    j = 2;
    a = R.traj{find(strcmp(allc, 'ES-QPID')), j};  b = R.traj{find(strcmp(allc, 'ES-noweight')), j};
    subplot(2,1,1);  hold on;
    patch([ton toff toff ton], [-0.05 -0.05 1.3 1.3], cgust, 'EdgeColor', 'none');
    h1 = plot(a.t, a.p, '-', 'Color', [0.10 0.10 0.10], 'LineWidth', 1.6);
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.1]);
    ylabel('P_{reg}');  title('(a) probability weight of the regulation state, L = 0.8 s');
    subplot(2,1,2);  hold on;
    patch([ton toff toff ton], [-0.1 -0.1 1.3 1.3], cgust, 'EdgeColor', 'none');
    plot([0 24], [1 1], 'k:', 'LineWidth', 1.0);
    h2 = plot(b.t, b.y, '--', 'Color', sty.ES_noweight{1}, 'LineWidth', 1.4);
    h3 = plot(a.t, a.y, '-',  'Color', sty.ES_QPID{1}, 'LineWidth', 2.0);
    box on;  grid on;  xlim([0 24]);  ylim([-0.05 1.3]);
    xlabel('Time (s)');  ylabel('Azimuth (deg)');  title('(b) response with and without the probability weights');
    legend([h2 h3], {'constant weight 0.5 (ES-noweight)', 'two-qubit weights (ES-QPID)'}, 'Location', 'southeast');
    savefig(fh, odir, 'fig_weights', dpi);
  end
  printf('Figures saved in %s\n', odir);
end

function fh = newfig(w, h)
  fh = figure('Visible', 'off', 'PaperUnits', 'inches', 'PaperPosition', [0 0 w h], ...
              'PaperSize', [w h], 'Position', [50 50 round(w*100) round(h*100)]);
  set(fh, 'DefaultAxesFontSize', 9, 'DefaultTextFontSize', 9, 'DefaultAxesFontName', 'Helvetica');
end

function savefig(fh, odir, name, dpi)
  print(fh, fullfile(odir, [name '.png']), '-dpng', sprintf('-r%d', dpi));
  close(fh);
end
