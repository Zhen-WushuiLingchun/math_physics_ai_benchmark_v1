import json, numpy as np
d = json.load(open('work/usage_messages.json', encoding='utf-8'))
msgs = [m for m in d['per_message'] if m['cost'] > 0]
A = np.array([[m['input'], m['cache_read'], m['output'], m['reasoning']] for m in msgs], float)
y = np.array([m['cost'] for m in msgs], float)
# least squares for 4 rates
coef, res, rank, sv = np.linalg.lstsq(A, y, rcond=None)
print('rows:', len(msgs))
print('rates per token: input=%.6e cache_read=%.6e output=%.6e reasoning=%.6e' % tuple(coef))
print('per million:     input=%.6f cache_read=%.6f output=%.6f reasoning=%.6f' % tuple(coef*1e6))
pred = A @ coef
print('max abs residual:', np.max(np.abs(pred - y)), ' rel:', np.max(np.abs(pred - y)/np.maximum(y,1e-12)))
# alternative: reasoning same price as output
A2 = np.array([[m['input'], m['cache_read'], m['output']+m['reasoning']] for m in msgs], float)
coef2, *_ = np.linalg.lstsq(A2, y, rcond=None)
print('alt(3-rate, reason=output): input=%.6f cache_read=%.6f output=%.6f  $/M' % tuple(coef2*1e6))
pred2 = A2 @ coef2
print('alt max abs residual:', np.max(np.abs(pred2 - y)))
