fileID = fopen('Gl_prm.txt');
formatSpec = '%s v1=%f v2=%f v3=%f %*s';
C=textscan(fileID,formatSpec);