from tests import OSmiPyTestCase

from osmipy import lexer, tokens
from osmipy.tokens import Token as Tok


class LexerTestCase(OSmiPyTestCase):

    def test_lexer(self):
        """Test that the lexer gives the good tokens"""

        tests = [
            ('cccccc', [Tok(tokens.ATOM, 'c')] * 6),
            ('C1=CC=CC=C1', [
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.DIGIT, 1),
                Tok(tokens.BOND, '='),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.BOND, '='),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.BOND, '='),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.DIGIT, 1),
            ]),
            ('N[C@](Br)(O)C', [
                Tok(tokens.ATOM, 'N'),
                Tok(tokens.LSPAR, '['),
                Tok(tokens.ATOM, 'C'),
                Tok(tokens.AT, '@'),
                Tok(tokens.RSPAR, ']'),
                Tok(tokens.LPAR, '('),
                Tok(tokens.ATOM, 'Br'),
                Tok(tokens.RPAR, ')'),
                Tok(tokens.LPAR, '('),
                Tok(tokens.ATOM, 'O'),
                Tok(tokens.RPAR, ')'),
                Tok(tokens.ATOM, 'C'),
            ])
        ]

        for s, r in tests:
            smiles_lexer = lexer.Lexer(s)
            lexed = [i for i in smiles_lexer.tokenize()]
            self.assertEqual(len(lexed) - 1, len(r))
            for i, token in enumerate(r):
                self.assertEqual(token.type, lexed[i].type, msg='{} of {}'.format(i, s))
                self.assertEqual(token.value, lexed[i].value, msg='{} of {}'.format(i, s))

            self.assertEqual(lexed[-1].type, tokens.EOF)
