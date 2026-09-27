"""Serve only the viewer directory plus the generated aircraft STEP download."""
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
import argparse
ROOT=Path(__file__).resolve().parent
class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs):super().__init__(*args,directory=str(ROOT),**kwargs)
    def do_GET(self):
        if self.path.split('?')[0]=='/downloads/spearhead.step':
            file=ROOT.parent/'output/spearhead.step'
            if not file.exists():self.send_error(404,'Generate the aircraft STEP first');return
            self.send_response(200);self.send_header('Content-Type','application/step');self.send_header('Content-Disposition','attachment; filename="spearhead.step"');self.send_header('Content-Length',str(file.stat().st_size));self.end_headers()
            with file.open('rb') as stream:self.copyfile(stream,self.wfile)
            return
        super().do_GET()
    def end_headers(self):
        self.send_header('Cache-Control','no-cache');self.send_header('X-Content-Type-Options','nosniff');super().end_headers()
    def list_directory(self,path):self.send_error(404);return None
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=4195);args=parser.parse_args()
    ThreadingHTTPServer(('0.0.0.0',args.port),Handler).serve_forever()
