classdef LabviewClient < handle
    
    properties
        socket
        server
        output
        visionSocket
        visionOutput
        visionInput
    end
    
    properties (Access = private)
        readTimeout
        data
    end
    
    methods
        
        function obj = LabviewClient()
            import java.io.*;
            import java.net.*;
            import java.util.ArrayList;
            import java.util.Collections;
            import java.util.List;
            import edu.ucsc.neurobiology.vision.util.*;
            
            % Create a server socket.
            disp('Creating server socket...');
            obj.server = ServerSocket(7887);
            
            host = InetAddress.getLocalHost();
            disp(['Serving on host: ', char(host.getHostName()), ' and port: 7887']);
            
            
            obj.socket = java.net.Socket();
            
            obj.readTimeout = 10;
        end
        
        function delete(obj)
            obj.close();
        end
        
        function connect(obj, port)
            import java.io.*;
            import java.net.*;
            import java.util.ArrayList;
            import java.util.Collections;
            import java.util.List;
            
            host = InetAddress.getLocalHost();
            obj.socket = Socket(host, port);
            obj.output = obj.socket.getOutputStream();
            
            % Send the spike finding command.
            obj.startSpikeFinding();
            
            % Listen for connections from Vision.
            obj.visionSocket = obj.server.accept();
            obj.visionOutput = obj.visionSocket.getOutputStream();
            obj.visionInput = obj.visionSocket.getInputStream();
        end
        
        function startSpikeFinding(obj)
            % Send the spike finding command.
            obj.output.write( java.nio.ByteBuffer.allocate(4).putInt(34).array() );
        end
        
        function initHandshake(obj)
            obj.visionOutput.write(java.nio.ByteBuffer.allocate(4).putInt(34).array());
            
            disp(obj.visionInput.read());
        end
        
        function readData(obj)
            fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data011/data011000.bin';
            fileID = fopen(fileName);
            % Read the whole file.
            obj.data = fread(fileID);
            fclose(fileID);
        end
        
%         function connect(obj, host, port)
%             % Connects to the specified host ip on the specified port.
%             
%             addr = java.net.InetSocketAddress(host, port);
%             timeout = 10000;
%             
%             try
%                 obj.socket.connect(addr, timeout);
%             catch x
%                 error(char(x.ExceptionObject.getMessage()));
%             end
%         end
        
        function streamData(obj)
            import java.io.*;
            import java.net.*;
            import java.util.ArrayList;
            import java.util.Collections;
            import java.util.List;
            import edu.ucsc.neurobiology.vision.util.*;
            
            for jj = 1 : length(obj.data)
                obj.visionOutput.write(obj.data(jj));
            end

%             obj.output.close();
%             obj.socket.close();
        end
        
        function close(obj)
            clear obj.data;
            if ~isempty(obj.output)
                obj.output.close();
            end
            obj.server.close();
            if ~isempty(obj.socket)
                obj.socket.close();
            end
            if ~isempty(obj.visionOutput)
                obj.visionOutput.close();
            end
            if ~isempty(obj.visionSocket)
                obj.visionSocket.close();
            end
            if ~isempty(obj.visionInput)
                obj.visionInput.close();
            end
        end
        
        function writeString(obj, v)
            try
                writer = java.io.PrintWriter(obj.socket.getOutputStream());
            catch x
                if isa(x, 'matlab.exception.JavaException')
                    error(char(x.ExceptionObject.getMessage()));
                end
                rethrow(x);
            end
            
            try
                writer.println(v);
%                 writer.write(v);
%                 writer.flush();
%                 writer.close();
            catch x
                if isa(x, 'matlab.exception.JavaException')
                    error(char(x.ExceptionObject.getMessage()));
                end
                rethrow(x);
            end
        end
        
        function write(obj, varargin)
            try
                stream = java.io.ObjectOutputStream(obj.socket.getOutputStream());
            catch x
                if isa(x, 'matlab.exception.JavaException')
                    error(char(x.ExceptionObject.getMessage()));
                end
                rethrow(x);
            end
            
            bytes = getByteStreamFromArray(varargin);
            
            try
                stream.writeObject(bytes);
            catch x
                if isa(x, 'matlab.exception.JavaException')
                    error(char(x.ExceptionObject.getMessage()));
                end
                rethrow(x);
            end
        end
        
        function fname = getFileName(obj, timeout)
            try
                fname = obj.listenToServer(timeout);
            catch
            end
            fname = char(fname);
        end
        
        function result = listenToServer(obj, timeout)
            % Timeout in seconds.
            in = obj.socket.getInputStream();
            
            start = tic;
            while in.available() == 0
                if timeout > 0 && toc(start) >= timeout
                    error('TcpConnection:ReadTimeout', 'Read timeout');
                end
            end
            
            stream = java.io.ObjectInputStream(in);
            
            result = stream.readObject();
            
%             varargout = getArrayFromByteStream(typecast(result, 'uint8'));
        end
        
        function result = read(obj)
            in = obj.socket.getInputStream();
            
            start = tic;
            while in.available() == 0
                if obj.readTimeout > 0 && toc(start) >= obj.readTimeout
                    error('TcpConnection:ReadTimeout', 'Read timeout');
                end
            end
            
            stream = java.io.ObjectInputStream(in);
            
            result = stream.readObject();
            
%             varargout = getArrayFromByteStream(typecast(result, 'uint8'));
        end
    end
end