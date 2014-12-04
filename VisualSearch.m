function VisualSearch(varargin)

global DIMS wRect XCENTER YCENTER STIM COLORS PICS VST w

% prompt={'SUBJECT ID' 'Session (1, 2, or 3)' 'Practice? 0 or 1'};
% defAns={'4444' '' ''};
% 
% answer=inputdlg(prompt,'Please input subject info',1,defAns);
% 
% ID=str2double(answer{1});
% % COND = str2double(answer{2});
% SESS = str2double(answer{2});
% prac = str2double(answer{3});


ID = 4444;
SESS = 1;

rng(ID); %Seed random number generator with subject ID
d = clock;

COLORS = struct;
COLORS.WHITE = [255 255 255];
COLORS.RED = [255 0 0];
COLORS.GREEN = [0 255 0];


DIMS = struct;
DIMS.grid_row = 4; %These have to be even numbers...
DIMS.grid_col = 4; %These have to be even numbers...
DIMS.grid_totes = DIMS.grid_row*DIMS.grid_col;
DIMS.maxside = 300;
DIMS.minside = 60;

STIM = struct;
STIM.trialdur = 3;
STIM.trials = 30;
STIM.blocks = 4;
STIM.totes = STIM.trials * STIM.blocks;

%% Find & load in pics
%find the image directory by figuring out where the .m is kept

[imgdir,~,~] = fileparts(which('VisualSearch.m'));

try
    cd([imgdir filesep 'IMAGES'])
catch
    error('Could not find and/or open the IMAGES folder.');
end

PICS =struct;
% if COND == 1;                   %Condtion = 1 is food. 
    PICS.in.lo = dir('good*.jpg');
    PICS.in.hi = dir('*bad*.jpg');
%     PICS.in.neut = dir('*water*.jpg');
% elseif COND == 2;               %Condition = 2 is not food (birds/flowers)
%     PICS.in.hi = dir('*bird*.jpg');
%     PICS.in.hi = dir('*flowers*.jpg');
%     PICS.in.neut = dir('*mam*.jpg');
% end
% picsfields = fieldnames(PICS.in);

picdiff = STIM.totes - length(PICS.in.lo);
if picdiff > 0;
    %if pic list is too short, increase it arbitrarily
    PICS.in.lo = [PICS.in.lo; PICS.in.lo(1:picdiff)];
elseif picdiff < 0;
    %Pic list is too long; arbitrarily increase it
    PICS.in.lo = PICS.in.lo(1:STIM.totes);
end

%Check if pictures are present. If not, throw error.
%Could be updated to search computer to look for pics...
if isempty(PICS.in.lo) || isempty(PICS.in.hi) %|| isempty(PICS.in.neut)
    error('Could not find pics. Please ensure pictures are found in a folder names IMAGES within the folder containing the .m task file.');
end

%% Setup Trial Variables
VST = struct;
VST.var.picnum_lo = reshape((randperm(length(PICS.in.lo))),STIM.trials,STIM.blocks);    %Pic random order of 120 hi cal images
VST.var.lo_loc = randi(DIMS.grid_totes,STIM.trials,STIM.blocks);                        %Location for lo-cal food.
for xx = 1:STIM.totes;
    %Load random sets of hi-cal food to fill in image grid. Each row in
    %VST.var.picnum_hi represents the 15 (or whatever #) of hi-cal foods to
    %display for each trial
    VST.var.picnum_hi(xx,1:(DIMS.grid_totes)) = randi(length(PICS.in.hi),1,16);
    VST.var.picnum_hi(xx,VST.var.lo_loc(xx)) = 0;
end
VST.data.rt = zeros(STIM.trials,STIM.blocks);
VST.data.correct = repmat(-999,STIM.trials,STIM.blocks);
VST.data.avg_rt = zeros(STIM.blocks,1);
VST.data.info.ID = ID;
% VST.data.info.cond = COND;               %Condtion 1 = Food; Condition 2 = animals
VST.data.info.session = SESS;
VST.data.info.date = sprintf('%s %2.0f:%02.0f',date,d(4),d(5));

%%
%change this to 0 to fill whole screen
DEBUG=0   ;

%set up the screen and dimensions

%list all the screens, then just pick the last one in the list (if you have
%only 1 monitor, then it just chooses that one)
Screen('Preference', 'SkipSyncTests', 1)

screenNumber= max(Screen('Screens'));

if DEBUG==1;
    %create a rect for the screen
    winRect=[0 0 640 480];
    %establish the center points
    XCENTER=320;
    YCENTER=240;
else
    %change screen resolution
    %Screen('Resolution',0,1024,768,[],32);
    
    %this gives the x and y dimensions of our screen, in pixels.
    [swidth, sheight] = Screen('WindowSize', screenNumber);
    XCENTER=fix(swidth/2);
    YCENTER=fix(sheight/2);
    %when you leave winRect blank, it just fills the whole screen
    winRect=[];
end

%open a window on that monitor. 32 refers to 32 bit color depth (millions of
%colors), winRect will either be a 1024x768 box, or the whole screen. The
%function returns a window "w", and a rect that represents the whole
%screen. 
[w, wRect]=Screen('OpenWindow', screenNumber, 0,winRect,32,2);

%%
%you can set the font sizes and styles here
Screen('TextFont', w, 'Arial');
%Screen('TextStyle', w, 1);
Screen('TextSize',w,30);

KbName('UnifyKeyNames');

%% Rects & other constants based off rects
%Set up rects for images with "DrawRectsGrid"
rects = DrawRectsGrid();
%Pre-set color frame size/location as rect + 10 for when an image is selected
% color_rects = [rects(1:2,:)-10; rects(3:4,:)+10];

%% Practice
% if prac == 1;
    DrawFormattedText(w,' Let''s practice.\n\nPress any key to continue.','center','center',COLORS.WHITE);
    Screen('Flip',w);
    KbWait([],2);
    
    DrawFormattedText(w,'In this task, you will see a grid of 16 images of food. It is your job to click on the image of the low calorie food as quickly as you can.\n\n Press any key to try a round.','center','center',COLORS.WHITE,60,[],[],1.5);
    Screen('Flip',w);
    KbWait([],2);
    
    while 1
        [prac_pics] = DrawPics4Trial(0,0);
        [~, prac_correct] = DoVisualSearch(0,0,prac_pics,rects);
        if prac_correct == 0;
            DrawFormattedText(w,'Try again!','center','center',COLORS.RED);
            Screen('Flip',w);
            WaitSecs(2);
        elseif prac_correct == 1;
            DrawFormattedText(w,'Very good!','center','center',COLORS.RED);
            Screen('Flip',w);
            WaitSecs(2);
            break
        end
    end
    
    DrawFormattedText(w,'In the real trials, you will only have 3 seconds to find the low calorie food,so you must move quickly! If you don''t find it in time, you will see "Time Expired" on the screen.\n\n Press any key when you are ready to move on to the task.','center','center',COLORS.WHITE,60,[],[],1.5);
    Screen('Flip',w);
    KbWait([],2);

    
% end


%% The Task

for block = 1:STIM.blocks;
    ibt = sprintf('Prepare for Block %d. \n\n\nPress any key when you are ready to begin.',block);
    DrawFormattedText(w,ibt,'center','center',COLORS.WHITE);
    Screen('Flip',w);
    KbWait();
    
    for trial = 1:STIM.trials;
        [pics] = DrawPics4Trial(trial,block);
        [VST.data.rt(trial,block), VST.data.correct(trial,block)] = DoVisualSearch(trial, block, pics, rects);
        %After trial, clear screen and wait in darkness for 500 ms.
        Screen('Flip',w);
        WaitSecs(.5);
    end
    
    %Interblock stats.
    Screen('Flip',w);   %clear screen first.
    
    block_text = sprintf('Block %d Results',block);
    %Find correct trials
    c = find(VST.data.correct(:,block) == 1);
    corr_count = sprintf('Number Correct:\t%d of %d',length(c),STIM.trials);  %Number correct = length of find(c)
    corr_per = length(c)*100/STIM.trials;                           %Percent correct = length find(c) / total trials
    corr_pert = sprintf('Percent Correct:\t%4.1f%%',corr_per);          %sprintf that data to string.
    
    if isempty(c)
        %Don't try to calculate avg RT, they got them all wrong (WTF?)
        %Display "N/A" for this block's RT.
        ibt_rt = sprintf('Average RT:\tUnable to calculate RT.');
    else
%         block_go = DPT.var.GoNoGo(:,block) == 1;                        %Find go trials
%         blockrts = DPT.data.rt(:,block);                                %Pull all RT data
        blockrts = VST.data.rt(c,block);                                         %Resample RT only if correct.
        avg_rt_block = fix(mean(blockrts)*1000);                        %Display avg rt in milliseconds.
        ibt_rt = sprintf('Average RT:\t\t\t%3d milliseconds',avg_rt_block);
    end
    
    ibt_xdim = wRect(3)/10;
    ibt_ydim = wRect(4)/4;

    DrawFormattedText(w,block_text,'center',wRect(4)/10,COLORS.WHITE);   %Next lines display all the data.
    DrawFormattedText(w,corr_count,ibt_xdim,ibt_ydim,COLORS.WHITE);
    DrawFormattedText(w,corr_pert,ibt_xdim,ibt_ydim+30,COLORS.WHITE);    
    DrawFormattedText(w,ibt_rt,ibt_xdim,ibt_ydim+60,COLORS.WHITE);
    %Screen('Flip',w);
    
    if block > 1
        % Also display rest of block data summary
        tot_trial = block * STIM.trials;
        totes_c = find(VST.data.correct==1);
        corr_count_totes = sprintf('Number Correct: \t%d of %d',length(totes_c),tot_trial);
        corr_per_totes = length(totes_c)*100/tot_trial;
        corr_pert_totes = sprintf('Percent Correct:\t%4.1f%%',corr_per_totes);
        
        if isempty(totes_c(totes_c ==1))
            %Don't try to calculate RT, they have missed EVERY SINGLE GO
            %TRIAL! 
            %Stop task & alert experimenter?
            tot_rt = sprintf('Block %d Average RT:\tUnable to calculate RT.',block);
        else
%             tot_go = DPT.var.GoNoGo == 1;
%             totrts = DPT.data.rt;
%             totrts = totrts(totes_c);
            tote_rts = VST.data.rt(c);
            avg_rt_tote = fix(mean(tote_rts)*1000);     %Display in units of milliseconds.
            tot_rt = sprintf('Average RT:\t\t\t%3d milliseconds',avg_rt_tote);
        end
        
        DrawFormattedText(w,'Total Results','center',ibt_ydim+120,COLORS.WHITE);
        DrawFormattedText(w,corr_count_totes,ibt_xdim,ibt_ydim+150,COLORS.WHITE);
        DrawFormattedText(w,corr_pert_totes,ibt_xdim,ibt_ydim+180,COLORS.WHITE);
        DrawFormattedText(w,tot_rt,ibt_xdim,ibt_ydim+210,COLORS.WHITE);
    end        
end

%% Save

% XXX: Need save function here.
end

%%
function [ rects ] = DrawRectsGrid(varargin)
%DrawRectGrid:  Builds a grid of squares with gaps in between.

global DIMS wRect XCENTER YCENTER

%Size of image will depend on screen size. First, an area approximately 80%
%of screen is determined. Then, images are 1/4th the side of that square
%(minus the 3 x the gap between images.

ylen = wRect(4)*8/10;           %Make square play area covering about 80% of vertical dimension of screen.
gap = 10;                       %Gap size between each image
square_side = fix((ylen - 3*gap)/4); %Size of image depends on size of screen.

squart_x = XCENTER-(ylen/2);
squart_y = YCENTER-(ylen/2);

rects = zeros(4,(DIMS.grid_totes));

for row = 1:DIMS.grid_row;
    for col = 1:DIMS.grid_col;
        currr = ((row-1)*DIMS.grid_col)+col;
        rects(1,currr)= squart_x + (row-1)*(square_side+gap);
        rects(2,currr)= squart_y + (col-1)*(square_side+gap);
        rects(3,currr)= squart_x + (row-1)*(square_side+gap)+square_side;
        rects(4,currr)= squart_y + (col-1)*(square_side+gap)+square_side;
    end
end

end

%%
function [trial_pics] = DrawPics4Trial(trial, block, varargin)
%Draw images for trial of Visual Search task;

global PICS VST w DIMS STIM

trial_pics = zeros(DIMS.grid_totes,1);

if trial == 0 && block == 0;
    %PRACTICE TRIAL.
    %Just pick first 15 hi cal & 1 st local food in list.
    for p = 1:16;
        thispic = imread(getfield(PICS,'in','hi',{p},'name'));
        trial_pics(p) = Screen('MakeTexture',w,thispic);
    end
    %Put lo cal food in square 10
    thispic = imread(getfield(PICS,'in','lo',{1},'name'));
    trial_pics(10) = Screen('MakeTexture',w,thispic);

else
    currtri = ((block-1)*STIM.trials)+trial;
    
%     trial_pics = zeros(DIMS.grid_totes,1);
    trial_picnums = VST.var.picnum_hi(currtri,:);
    
    
    for p = 1:DIMS.grid_totes;
        if trial_picnums(p) == 0;
            %load lo-cal food
            thispic = imread(getfield(PICS,'in','lo',{VST.var.picnum_lo(trial,block)},'name'));
        else
            %load hi-cal food
            thispic = imread(getfield(PICS,'in','hi',{trial_picnums(p)},'name'));
            
        end
        trial_pics(p) = Screen('MakeTexture',w,thispic);
    end
end


end

%%
function [rt, correct] = DoVisualSearch(trial, block, pics, rects, varargin)

global VST w COLORS DIMS STIM

correct = -999;

Screen('DrawTextures',w,pics,[],rects);
startRT = Screen('Flip',w);
telap = GetSecs - startRT;

if trial == 0 && block == 0;
    trial_duration = 5000;     %If practice, duration = some arbitrary large time.
    boxcheck = 10;             %If practice, lo cal food = box 10.
else
    trial_duration = STIM.trialdur;
    boxcheck = VST.var.lo_loc(trial,block);
end

while telap < trial_duration;
    telap = GetSecs - startRT;
    [x,y,button] = GetMouse();
    if button(1);
        %test if mouse clicked in one of our precious boxes.
        %This creates arrays of 1s or 0s for > or < min/max dimensions
        rt = GetSecs - startRT;
        
        xmin = rects(1,:)<=x;
        xmax = rects(3,:)>=x;
        ymin = rects(2,:)<=y;
        ymax = rects(4,:)>=y;
        
        %This tests if there are any cases where all of the above are true.
        clickedonbox = find(xmin & xmax & ymin & ymax);
        
        if ~isempty(clickedonbox);
            rect_zoom = rects(:,clickedonbox);
            %                 color_zoom = color_rects(:,clickedonbox);
            side = 0;
            
            %Which quadrant is it in?
            switch clickedonbox;
                case {1,2,5,6}
                    quad = 1;
                case {3,4,7,8}
                    quad = 2;
                case{9,10,13,14}
                    quad = 3;
                case{11,12,15,16}
                    quad = 4;
            end
            if clickedonbox == boxcheck;
                %They have clicked lo-cal food. Do the zoom thing.
                %Depending on image location, will need to pin image from
                % different corners of box.
                % Note: These are set up for a 4x4 grid and are not
                % adaptable to other arrangements (currently).
                correct = 1;
                break
            else
                %This person mis-clicked on a hi-cal food
                correct = 0;
                break
            end
        else
            FlushEvents();
        end
    end
end

if correct == -999
    %If correct = -999, then no press was recorded. Throw incorrect response.
    %XXX: Should there be some zoomy thing?
    DrawFormattedText(w,'Time Expired','center','center',COLORS.RED);
    rt = -999;
    correct = 0;
    
elseif correct == 1;
    while side < DIMS.maxside
        side = side + 10;
        quad_zoom = [rect_zoom+[0;0;side;side] rect_zoom+[0;-side;side;0] rect_zoom+[-side;0;0;side] rect_zoom+[-side;-side;0;0]];
        quadc_zoom = quad_zoom + repmat([-10;-10;10;10],1,4);
        Screen('DrawTextures',w,pics,[],rects);
        Screen('FillRect',w,COLORS.GREEN,quadc_zoom(:,quad));
        Screen('DrawTexture',w,pics(clickedonbox),[],quad_zoom(:,quad));
        Screen('Flip',w);
    end
    WaitSecs(.5);
elseif correct == 0;
    while side < DIMS.minside
        side = side + 1;
        quad_zoom = rect_zoom + [side;side;-side;-side];
        %         quad_zoom = [rect_zoom+[0;0;-side;-side] rect_zoom+[0;side;-side;0] rect_zoom+[side;0;0;-side] rect_zoom+[side;side;0;0]];
        quadc_zoom = rect_zoom + [-10;-10;10;10];
        Screen('DrawTextures',w,pics,[],rects);
        Screen('FillRect',w,COLORS.RED,quadc_zoom);
        Screen('DrawTexture',w,pics(clickedonbox),[],quad_zoom(:,quad));
        Screen('Flip',w);
    end
    WaitSecs(.5);
end


FlushEvents();    

end

