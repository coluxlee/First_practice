%% =====================================================================
%  HH 모델 과제 Part A — 발표용 코드 (HH_present.m)
%  [기반] HH 엔진 → A-1 → A-2 → A-3 → A-4
%  설명하기 쉽게 [기반](함수)을 맨 앞에 둔 발표용 배치
%  
%  단위 [ 시간 ms · 전압 mV · 전류 uA/cm^2 · 전도도 mS/cm^2 · 용량 uF/cm^2 ]
% =====================================================================

%% ■■■■■■■■■■■■■■■■■■■■■■ [기반] HH 엔진 ■■■■■■■■■■■■■■■■■■■■■■
% A-2 ~ A-4는 모두 이 엔진을 불러서 씀
%   기반 ① α·β 함수 6개 : 게이트가 열리고 닫히는 속도 (교재 12.43–12.46)
%   기반 ② run_hh       : 상수 + 조건 + 공책 + 오일러 반복 = 시뮬레이션 1회

%% 기반 ① 게이트 속도 함수 α·β 6개 (PPT 7쪽)
% 입력 : 지금 막전위 Vm [mV]  →  출력 : 속도 [1/ms]
% 교재 식의 V = '휴지전위로부터의 변위' (탈분극이면 V < 0) → 입구에서 V = V_rp - Vm 으로 변환
% 함수 = 밀폐된 방 : 바깥 상수가 안 보임 → V_rp를 방마다 다시 적음

function a = alpha_n(Vm)                
    V_rp = -60;  V = V_rp - Vm;         
    if abs(V + 10) < 1e-6               
        a = 0.1;                        
    else
        a = 0.01 * (V + 10) / (exp((V + 10)/10) - 1);
    end
end

function b = beta_n(Vm)                 
    V_rp = -60;  V = V_rp - Vm;
    b = 0.125 * exp(V/80);
end

function a = alpha_m(Vm)                
    V_rp = -60;  V = V_rp - Vm;
    if abs(V + 25) < 1e-6               
        a = 1.0;                        
    else
        a = 0.1 * (V + 25) / (exp((V + 25)/10) - 1);
    end
end

function b = beta_m(Vm)                 
    V_rp = -60;  V = V_rp - Vm;
    b = 4 * exp(V/18);
end

function a = alpha_h(Vm)                
    V_rp = -60;  V = V_rp - Vm;
    a = 0.07 * exp(V/20);
end

function b = beta_h(Vm)                 
    V_rp = -60;  V = V_rp - Vm;
    b = 1 / (exp((V + 30)/10) + 1);     
end

%% 기반 ② HH 엔진 run_hh : 시뮬레이션 1회
% 입구 : I_amp = 펄스 세기 [uA/cm^2],  pw = 펄스 폭 [ms]
% 출구 : Vm, t, gNa, gK, Iin (모두 칸 N개짜리 공책 → 그래프 재료)
% 부를 때마다 빈 방에서 새로 시작 → 몇 번을 불러도 결과가 섞이지 않음

function [Vm, t, gNa, gK, Iin] = run_hh(I_amp, pw)  % 전류/펄스폭(ms)

    % [1] 상수
    gNa_bar = 120;  gK_bar = 36;  gL = 0.3;     % 최대 전도도 [mS/cm^2] (_bar = 문이 다 열렸을 때)
    E_Na = 55;  E_K = -72;  E_L = -49.4;        % 평형전위 [mV]
    C_m = 1;  V_rp = -60;                       % 막 용량 [uF/cm^2], 휴지전위 [mV]

    % [2] 조건1 시간축
    dt = 0.01;              
    T  = 30;                
    t  = 0:dt:T;          
    N  = length(t);        

    % [3] 조건2 입력 전류
    t_on  = 5;                                  
    k_on  = round(t_on/dt) + 1;                 
    k_off = min(k_on + round(pw/dt), N);      
    Iin   = zeros(1, N);                 
    Iin(k_on : k_off-1) = I_amp;     

    % [4] 조건3 게이트변수 초기값
    Vm  = zeros(1, N);  m  = zeros(1, N);  h = zeros(1, N);  n = zeros(1, N);  
    gNa = zeros(1, N);  gK = zeros(1, N);                                     
    Vm(1) = V_rp;                                                  
    m(1) = alpha_m(Vm(1)) / (alpha_m(Vm(1)) + beta_m(Vm(1)));     
    h(1) = alpha_h(Vm(1)) / (alpha_h(Vm(1)) + beta_h(Vm(1)));    
    n(1) = alpha_n(Vm(1)) / (alpha_n(Vm(1)) + beta_n(Vm(1)));    

    % [5] 오일러 반복 : k번째 칸 → k+1번째 칸 
    for k = 1:N-1                             

        gNa(k) = gNa_bar * m(k)^3 * h(k);
        gK(k)  = gK_bar  * n(k)^4;

        I_Na = gNa(k) * (Vm(k) - E_Na);            
        I_K  = gK(k)  * (Vm(k) - E_K);
        I_L  = gL     * (Vm(k) - E_L);
        % 기울기
        dVdt = (Iin(k) - I_Na - I_K - I_L) / C_m;  
        dmdt = alpha_m(Vm(k)) * (1 - m(k)) - beta_m(Vm(k)) * m(k); 
        dhdt = alpha_h(Vm(k)) * (1 - h(k)) - beta_h(Vm(k)) * h(k); 
        dndt = alpha_n(Vm(k)) * (1 - n(k)) - beta_n(Vm(k)) * n(k);
        % 오일러 갱신
        Vm(k+1) = Vm(k) + dVdt * dt;
        m(k+1)  = m(k)  + dmdt * dt;
        h(k+1)  = h(k)  + dhdt * dt;
        n(k+1)  = n(k)  + dndt * dt;
    end

    % 마지막 칸 : for문이 N-1까지라 전도도 N번 칸이 비어 있음 → 직접 채움
    gNa(N) = gNa_bar * m(N)^3 * h(N);
    gK(N)  = gK_bar  * n(N)^4;
end

%% ■■■ A-1 : 유도한 식은 엔진의 어디에? ■■■
% 키르히호프 전류 법칙 : I_in = I_C + I_Na + I_K + I_L,  I_C = C_m dVm/dt
%   → C_m dVm/dt = I_in - I_Na - I_K - I_L
% 코드 : 기반 ② [5] (c)의 한 줄 →  dVdt = (Iin(k) - I_Na - I_K - I_L) / C_m;
% (손으로 유도하는 문제라 추가 코드 없음 — 엔진 안에 이미 들어 있음)

%% ■■■■■■■■■■■■■■■■■■■■■■ 실행부 : A-2 ~ A-4 ■■■■■■■■■■■■■■■■■■■■■■
clear; clc; close all;          % 이전 실행의 변수 · 명령창 · 그림 창 정리

%% ■■■ A-2 : 1 uA vs 20 uA (펄스 폭 1 ms) ■■■
% 추가한 것 : 엔진을 세기만 바꿔 두 번 부르기 + 그림 1장(위아래 2칸) + 숫자 출력
% 펄스 모양은 기반 ② [3]이 만듦 → 여기선 run_hh(세기, 폭)의 값만 바꿔서 부름
% 출구 5개는 '순서대로' 바깥 이름에 담김 → 이름 끝 1, 20 = "몇 uA 실험의 공책인지" 표시
pw = 1;                                          % 펄스 폭 [ms] (함수 안의 pw와는 다른 상자)
[Vm1,  t, gNa1,  gK1,  Iin1 ] = run_hh(1,  pw);  
[Vm20, t, gNa20, gK20, Iin20] = run_hh(20, pw); 

% 그림 1 : 칸마다 왼쪽 축(자동 파랑) = 결과(막전위), 오른쪽 축(자동 주황) = 원인(넣어 준 전류)
% 두 칸의 축 범위를 똑같이 고정(ylim) → 위아래를 눈으로 비교 가능 / LineWidth = 선 굵기
figure(1);                                       % 1번 그림 창
subplot(2,1,1);                                  % 2행 1열로 나눈 칸 중 1번(위) : 1 uA
yyaxis left;  plot(t, Vm1,  'LineWidth', 1.5);  ylabel('V_m [mV]');      ylim([-80 60]);
yyaxis right; plot(t, Iin1, 'LineWidth', 1.0);  ylabel('I_{in} [\muA]'); ylim([0 25]);
xlabel('Time [ms]'); title('A-2 : I_{in} = 1 \muA, 1 ms'); grid on;

subplot(2,1,2);                                  % 2번(아래) 칸 : 20 uA (구조는 위와 같음)
yyaxis left;  plot(t, Vm20,  'LineWidth', 1.5); ylabel('V_m [mV]');      ylim([-80 60]);
yyaxis right; plot(t, Iin20, 'LineWidth', 1.0); ylabel('I_{in} [\muA]'); ylim([0 25]);
xlabel('Time [ms]'); title('A-2 : I_{in} = 20 \muA, 1 ms'); grid on;

% 숫자 출력 : %.1f = 소수 첫째 자리 숫자가 들어갈 빈칸 (뒤의 값이 순서대로 채움), \n = 줄바꿈
fprintf('A-2  1 uA : Vm 최대 = %.1f mV\n', max(Vm1));
fprintf('A-2 20 uA : Vm 최대 = %.1f mV, 최소 = %.1f mV\n', max(Vm20), min(Vm20));

%% ■■■ A-3 : 같은 조건(1 uA · 20 uA)의 G_Na, G_K ■■■
% 추가한 것 : 그림 1장(위아래 2칸) + 숫자 출력 — 계산 코드는 0줄
% 엔진이 매 칸 gNa(k), gK(k)를 공책에 적어 둠 → 출구로 받은 공책을 그리기만 함
% (칸 번호 없이 gNa = ... 로 쓰면 매번 덮어써서 마지막 값 1개만 남음 → 그래프 불가)
% plot(가로, 세로, 색) 세트 2개 = 선 2개 : gNa 자홍('m'), gK 초록('g') / legend는 그린 순서대로 이름 붙음
figure(2);                                       % 2번 그림 창
subplot(2,1,1);                                  % 위 칸 : 1 uA (대조군) ※ 세로축 0~0.4 = 아래 칸의 약 1/100
plot(t, gNa1, 'm', t, gK1, 'g', 'LineWidth', 1.5);
ylabel('Conductance [mS/cm^2]'); legend('G_{Na}', 'G_K');
xlabel('Time [ms]'); title('A-3 : I_{in} = 1 \muA'); grid on;

subplot(2,1,2);                                  % 아래 칸 : 20 uA (Na+ 먼저 → K+ 나중)
plot(t, gNa20, 'm', t, gK20, 'g', 'LineWidth', 1.5);
ylabel('Conductance [mS/cm^2]'); legend('G_{Na}', 'G_K');
xlabel('Time [ms]'); title('A-3 : I_{in} = 20 \muA'); grid on;

fprintf('A-3 20 uA : G_Na 최대 = %.1f, G_K 최대 = %.1f mS/cm^2\n', max(gNa20), max(gK20));

%% ■■■ A-4 : Strength-duration curve (펄스 폭 0.1 ~ 2 ms) ■■■
% 추가한 것 : 이분탐색 (반복문 안의 반복문) + 그림 1장
% 목표 : 펄스 폭마다 AP를 '겨우' 일으키는 최소 전류(문턱) 찾기
% 발화 판정 : Vm이 한 번이라도 0 mV를 넘으면 AP (발화는 +27 mV 이상, 불발은 -47 mV 이하 → 0에서 뚜렷이 갈림)
widths = 0.1:0.1:2;                         % 시험할 펄스 폭 20개 [ms]
I_th   = zeros(1, length(widths));          % 찾은 문턱을 적을 공책 : widths(w) ↔ I_th(w) 한 쌍

for w = 1:length(widths)                    % 바깥 반복 : 펄스 폭 하나씩
    lo = 0;                                 % 확실히 불발하는 전류 (0 = 자극 없음)
    hi = 200;                               % 확실히 발화하는 전류
    if max(run_hh(hi, widths(w))) <= 0      % 안전장치 : hi에서도 불발이면 탐색 범위가 틀린 것
        warning('펄스 폭 %.1f ms에서 200 uA로도 불발 → hi를 더 크게 잡을 것', widths(w));
    end
    for iter = 1:20                         % 안쪽 반복 : 범위를 20번 반으로 (200/2^20 ≈ 0.0002 uA)
        mid = (lo + hi) / 2;                % 가운데 값을 시험
        Vm_test = run_hh(mid, widths(w));   % 출구를 1개만 받으면 첫 번째(Vm)만 들어옴
        if max(Vm_test) > 0                 % 발화 → 문턱은 mid 이하 → 위쪽 절반 버림
            hi = mid;
        else                                % 불발 → 문턱은 mid 초과 → 아래쪽 절반 버림
            lo = mid;
        end
    end
    I_th(w) = hi;                           % 발화하는 가장 작은 값 = 문턱 → w번 칸에 기록
    fprintf('A-4  펄스 폭 %.1f ms : 문턱 전류 = %.2f uA\n', widths(w), I_th(w));   % 진행 상황
end

figure(3)
plot(widths, I_th, 'o-', 'LineWidth', 1.5)  % 'o-' : 점(o)을 찍고 선(-)으로 이음
xlabel('Pulse width [ms]'); ylabel('Threshold current [\muA]')
title('A-4: Strength-duration curve'); grid on


 %깃허브용 추가 글
